struct tree_node {
	size_t down;		/* To first daughter node or target part */
	size_t next;		/* To next subnode of father or to unkle */
	float pos[3];		/* Center of node */
	long npart;		/* Number of particles in node */
	float size;		/* Spatial extent of node */
} *tree;

size_t treesize;		/* Number of nodes in tree */


void build_tree();
void init_tree();
int pos2idx(size_t, size_t);
