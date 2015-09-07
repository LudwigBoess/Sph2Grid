#include "allvars.h"
#include "tree.h"

#ifdef HSMLFIND

#define NNGB 35
#define DNNGB 3
#define NGBMAX 50		/* Length of NGBlist */
#define MAXHSMLITERATIONS 10
/*  */
#define KERNELCONST 2.546479089

/* 0.5*sqrt(1+1/sin(pi/4)^2) =  */
#define MAXCUBEDISTFAC 0.86602544784545

size_t ngblist[NGBMAX + 1] = { 0 }, ngbcnt = 0;

float iterate_hsml_bisect(size_t);
void ngbfind_open_node(size_t, size_t);
void find_ngb_tree(size_t, float, size_t);
float guess_hsml(size_t, size_t);
float rho_SPH(size_t, float);
float W(size_t, size_t, float);

void find_hsml()
{
	size_t ipart;
	double tmp;

	fprintf(stdout, "Starting Neighbour finding ... ");

	for (ipart = 0; ipart < ThisTask.Npart[1]; ipart++) {
		tmp = iterate_hsml_bisect(ipart);
		printf("ipart=%zu newhsml=%g gadgethsml=%g delta=%g \n",
		       ipart, tmp, Gas[ipart].Hsml,
		       (tmp - Gas[ipart].Hsml) / Gas[ipart].Hsml);
		Gas[ipart].Hsml = iterate_hsml_bisect(ipart);
		Gas[ipart].Rho = rho_SPH(ipart, Gas[ipart].Hsml);
	}

	printf("done\n");

	return;
}


/* Find hsml using bisection */
float iterate_hsml_bisect(size_t ipart)
{
	int cnt = 0;
	float hsml = 0, left = 0, right = 0;

	hsml = guess_hsml(ipart, 0);

	ngbcnt = 0;

	for (;;) {

		find_ngb_tree(ipart, hsml, 0);

		if (!((ngbcnt > NNGB + DNNGB) || (ngbcnt < NNGB - DNNGB)))
			break;

		if (ngbcnt < NNGB - DNNGB) {
			left = fmax(hsml, left);
		} else {
			if (right != 0) {
				if (hsml < right)
					right = hsml;
			} else {
				right = hsml;
			}
		}

		if (left == 0 && right > 0) 
			hsml /= 1.26;	/* half volume */

		if (left > 0 && right == 0) 
			hsml *= 1.26;	/* double volume */

		if (left > 0 && right > 0)	/* average volume */
			hsml = pow(0.5 * (p3(left) + p3(right)), 1. / 3);

		if (cnt > MAXHSMLITERATIONS) {
			printf("\n Max Number of Iterations reached <%zu>\n",
			     ipart);
			break;
		}

		cnt++;

	}

	return (hsml);
}

/* Find all neighbours closer than hsml via the tree*/
void find_ngb_tree(size_t ipart, float hsml, size_t node)
{
	int i;
	size_t son = 0, jpart = 0;
	float d = 0, nx, ny, nz;
#ifdef PERIODIC
	int k;
	double boxsize = Snap.Boxsize;
#endif

	if (!node)
		ngbcnt = 0;

	if (tree[node].npart > 1) {	/* Decline into proper nodes */
		for (i = 0; i < 8; i++) {
			son = tree[node].down + i;

#ifdef PERIODIC
			for (k = 0; k < 8; k++) {	/* Mirror Particle */
				nx = !(k & 0x1) ?
				    tree[son].pos[0] : (tree[son].pos[0] > 0 ? 
                            tree[son].pos[0] - boxsize :
							tree[son].pos[0] +boxsize);
				ny = !(k & 0x2) ? 
                    tree[son].pos[1] : (tree[son].pos[1] > 0 ? 
                        tree[son].pos[1] - boxsize : 
                        tree[son].pos[1] + boxsize);
				nz = !(k & 0x4) ? 
                    tree[son].pos[2] : (tree[son].pos[2] > 0 ? 
                        tree[son].pos[2] - boxsize : 
                        tree[son].pos[2] + boxsize);
#else
			nx = tree[son].pos[0];
			ny = tree[son].pos[1];
			nz = tree[son].pos[2];
#endif

			d = sqrt(p2(P[ipart].Pos.x - nx) 
				 + p2(P[ipart].Pos.y - ny)
				 + p2(P[ipart].Pos.z - nz);

			if (d - 0.5 * sqrt3 * tree[son].size - hsml < 0)
				find_ngb_tree(ipart, hsml, son);
#ifdef PERIODIC
		}
#endif
    }
    } else if (tree[node].npart == 1) {	/* leaf */
        jpart = tree[node].down;
    }

#ifdef PERIODIC
    for (k = 0; k < 8; k++) {	/* Mirror Particle */
	    nx = !(k & 0x1) ?
	        tree[son].pos[0] : (tree[son].pos[0] > 0 ?
			    	tree[son].pos[0] - boxsize : 
                    tree[son].pos[0] + boxsize);
	    ny = !(k & 0x2) ? 
            tree[son].pos[1] : (tree[son].pos[1] > 0 ? 
                    tree[son].pos[1] - boxsize : 
                    tree[son].pos[1] + boxsize);
	    nz = !(k & 0x4) ? 
            tree[son].pos[2] : (tree[son].pos[2] > 0 ? 
                    tree[son].pos[2] - boxsize : 
                    tree[son].pos[2] + boxsize);
#else
    nx = tree[son].pos[0];
    ny = tree[son].pos[1];
    nz = tree[son].pos[2];
#endif

    d = sqrt(p2(P[ipart].Pos.x - nx)
	 + p2(P[ipart].Pos.y - ny) 
	 + p2(P[ipart].Pos.z - nz));

    if (d < hsml && ngbcnt < NGBMAX) {
	    ngblist[ngbcnt] = jpart;	/* Add neighbour */
    	ngbcnt++;
    }
#ifdef PERIODIC
    }
#endif
}

return;
}

/* 
 * Find a first guess from the tree recursively 
 * by estimating the volume per particle from
 * the local number density and compute an hsml 
 * for NNGB neighbours 
 */
float guess_hsml(size_t ipart, size_t node)
{
	size_t son;
	float hsml;

	if (tree[node].npart >= NNGB) {
		son = tree[node].down + pos2idx(ipart, node);
		hsml = guess_hsml(ipart, son);
	} else {
		hsml = tree[node].size	/* MAGIC :-D */
		    * pow(NNGB / tree[node].npart / fourpithirds, 1. / 3.);
	}
	return (hsml);
}

/* Springel 2005 eq. 5*/
float rho_SPH(size_t ipart, float hsml)
{
	size_t j, jpart;
	float sum = 0;

	find_ngb_tree(ipart, hsml, 0);

	for (j = 1; j < ngbcnt; j++) {
		jpart = ngblist[j];
		sum += P[jpart].Mass * W(ipart, jpart, hsml);
	}

	return (sum);
}

/* 3D Kernel (Springel 05) */
float W(size_t i, size_t j, float h)
{
	float r, d, result;

	r = sqrt(pow(P[i].Pos.x - P[j].Pos.x, 2)
		 + pow(P[i].Pos.y - P[j].Pos.y, 2)
		 + pow(P[i].Pos.z - P[j].Pos.z, 2));

	d = r / h;

	if (d <= 0.5) {
		result = 1 - 6 * d * d + 6 * d * d * d;
	} else if (d < 1) {
		result = 2 * (1 - d) * (1 - d) * (1 - d);
	} else {
		return (0);
	}
	return (KERNELCONST / (h * h * h) * result);
}
#endif
