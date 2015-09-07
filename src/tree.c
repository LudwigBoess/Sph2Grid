#include "allvars.h"
#include "tree.h"

#define NODES_PER_PARTICLE 4.4

size_t refine();
int idx2sign(int, int);

size_t npart, max_tree_size;

void init_tree()
{
	int i;

	my_assert(ThisTask.NTask == 1., "Tree not parallelised yet");

	npart = 0;
	for (i = 0; i < N_part_types; i++)
		npart += ThisTask.Npart[i];

	tree = (struct tree_node *) malloc(NODES_PER_PARTICLE *
					   ThisTask.Npart[0] *
					   sizeof(struct tree_node));

	tree[0].down = -1;
	tree[0].pos[0] = 0;	/* Box is already centered */
	tree[0].pos[1] = 0;
	tree[0].pos[2] = 0;
	tree[0].npart = 0;
	tree[0].size = (float) Snap.Boxsize;

	treesize = 1;
	max_tree_size = NODES_PER_PARTICLE * ThisTask.Npart[0];

	return;
}

void build_tree()
{
	size_t ipart, node, idx;
	int part_done;


	fprintf(stdout, "Starting Tree construction \n");

	for (ipart = 0; ipart < npart; ipart++) {
		node = 0;
		part_done = 0;

		while (!part_done) {
			tree[node].npart++;
			switch (tree[node].npart) {
			case 1:	/* empty, sort in */
				tree[node].down = ipart;
				part_done++;
				break;
			case 2:	/* full, but leaf, refine */
				node = refine(node, ipart);
				break;
			default:	/* full and not leaf, decline */
				idx = pos2idx(ipart, node);
				node = tree[node].down + idx;
				break;
			}
		}
	}

	fprintf(stdout, "Tree construction done\n"
		"Using <%Zu> treenodes for <%Zu> particles \n", treesize,
		npart);


	return;
}

/* Refines parent node, puts old particle
 * in place and returns proper son for ipart*/
size_t refine(size_t parent, size_t ipart)
{
	size_t son, target, oldpart;
	int i;

	oldpart = tree[parent].down;	/* Remember old particle */

	tree[parent].down = treesize;

	treesize += 8;		/* Add 8 more nodes */

	my_assert(treesize < max_tree_size, "Max tree size reached");

	for (i = 0; i < 8; i++) {	/* Init new son nodes */
		son = tree[parent].down + i;

		tree[son].down = -1;
		tree[son].npart = 0;
		tree[son].size = tree[parent].size / 2.;

		tree[son].pos[0] = tree[parent].pos[0]
		    + idx2sign(i, 0) * 0.5 * tree[son].size;
		tree[son].pos[1] = tree[parent].pos[1]
		    + idx2sign(i, 1) * 0.5 * tree[son].size;
		tree[son].pos[2] = tree[parent].pos[2]
		    + idx2sign(i, 2) * 0.5 * tree[son].size;
	}

	/* Treat old particle */
	target = tree[parent].down + pos2idx(oldpart, parent);
	tree[target].down = oldpart;
	tree[target].npart++;

	/* Find new target */
	target = tree[parent].down + pos2idx(ipart, parent);

	return (target);	/* Return new target node for ipart */
}

/* Find index of subnode from relative position */
int pos2idx(size_t ipart, size_t node)
{
	float dx, dy, dz;

	dx = P[ipart].Pos[0] - tree[node].pos[0];
	dy = P[ipart].Pos[1] - tree[node].pos[1];
	dz = P[ipart].Pos[2] - tree[node].pos[2];

	return ((dx > 0) + ((dy > 0) << 1) + ((dz > 0) << 2));
}

/* Convert node index to sign of component comp */
int idx2sign(int idx, int comp)
{
	return (-1 + 2 * ((idx & (1 << comp)) >> comp));
}
