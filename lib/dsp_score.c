#!/usr/bin/env runc
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

typedef struct {
    char cluster[256];
    char gene[256];
    int domain;
    int begin;
    int end;
} DOMAIN;

/* Prototype Declaration */
int search_boundary_i(char **alignment, int n_alignment, char **domain, char **cluster_flgs, int n_clusters, int i, int search_from, int search_to);
int search_boundary_i_local(char **alignment, int n_alignment, char **domain, char **cluster_flgs, int n_clusters, int i, int search_from, int search_to, double *max_score);
int search_boundary_i_local_fine(char **alignment, int n_alignment, char **domain, char **cluster_flgs, int n_clusters, int i, int search_from, int search_to, double *max_score);
int search_boundary_i_local_step(char **alignment, int n_alignment, char **domain, char **cluster_flgs, int n_clusters, int i, int search_from, int search_to, int step, double *max_score);
int search_boundary_local (char **alignment, int n_alignment, char **domain, char **cluster_flgs, int n_clusters, int search_from, int search_to, double *max_score);
int search_boundary (char **alignment, int n_alignment, char **domain, char **cluster_flgs, int n_clusters, int search_from, int search_to, double *max_score);
double calc_score_for_line (char **alignment, int n_alignment, char **domain, char **cluster_flgs, int n_clusters, int i);
double calc_score_for_lines (char **alignment, int n_alignment, char **domain, char **cluster_flgs, int n_clusters);
double calc_score(char **alignment, int n_alignment, char **domain, char **cluster_flgs, int n_clusters) ;
double score_between (int n_alignment, char **alignment, char **domain, char **gap, int n_clusters, char **cluster_flgs);
double sum_of_pairs_for_line_local(char **a, char **d, char **cluster_flgs, int n, int n_clusters, int i, int from, int to);
double sum_of_pairs_between(char **a, char **d, char **g, char **cluster_flgs, int n, int n_clusters, int i, int i2);
int count_gaps(char **alignment, char **domain, char **cluster_flgs, int i, int n, int cluster);
int count_gaps_pair(char **alignment, char **domain, char **cluster_flgs, int i, int i2, int n, int cluster);
int count_gaps_local(char **alignment, char **domain, char **cluster_flgs, int i, int n, int from, int to);
int count_gaps_local_pair(char **alignment, char **domain, char **cluster_flgs, int i, int i2, int n, int from, int to);
int count_gaps_local_with_flg(char **alignment, char **domain, char **cluster_flgs, int i, int n, int from, int to);
int count_gaps_between(char **alignment, char **domain, char **cluster_flgs, int i, int i2, int n, int n_clusters, int cluster1, int cluster2);
int eval_gap(char **alignment, char **domain, char **cluster_flgs, int i, int n, int start, int end, int cluster);
int eval_gap_with_flg(char **alignment, char **domain, char **cluster_flgs, int i, int n, int start, int end, int cluster);
int find_aa(char alignment[], char domain[], int start, int end, int cluster);
void set_boundary(char **alignment, int n_alignment, char **domain, char **cluster_flgs, int search_from, int search_to, int j_boundary);
void set_boundary_i(char **alignment, char **domain, int i, int j_boundary, int search_from, int search_to);
void set_boundary_i_slow(char **alignment, char **domain, int i, int j_boundary, int search_from, int search_to);
void map_on_alignment(DOMAIN *a_domain, char **alignment, char **gene, int n_alignment, char **domain, char **cluster_flgs, int cluster);
void map_gaps (char *alignment, char *gap, int n_pos);
int count_aa (char **a, char **d, int n);
int count_aa_i (char **a, char **d, int i);
int get_cluster (DOMAIN *a_domain, char **clusters, int n_clusters);
void get_score_matrix(int score[][256], char *file);
int read_matrix (char *file, char key[], int mat[][256]);
int read_alignment (char *file, char **alignment, char **gene, int n_alignment);
int count_alignment (char *file);
int count_dclst (char *file);
int read_dclst (char *file, DOMAIN dclst[]);
void *my_calloc (size_t n, size_t size);
void *my_realloc (void *ptr, size_t size);
void *my_malloc (size_t size);
void chomp(char *s);
int split_int(char *line, int array[]);
int split_csv_to_int(char *line, int *array);
int split_char(char *line, char array[]);
int split_line (char *line, char f[][BUFSIZ]);
int split_line_to_domain(char *line, DOMAIN *dclst);
void print_char(char array[], int n);
void print_int(int array[], int n);
void print_mat(int mat[][256], int n);
int max (int x, int y);
int min (int x, int y);
/* Prototype Declaration End */

/* settings */
double GAP_EXT = 0.5;
double GAP_OPEN = 10;

/* score matrix */
int SCORE[256][256];

int main(int argc, char *argv[])
{
    char *program = argv[0];
    char usage[BUFSIZ];
    sprintf(usage,   "Usage: %s -m SCORE_MATRIX -t DCLST_FILE -a ALIGNMENT\n", program);
    sprintf(usage, "%s-e GAP_EXT\n", usage);
    sprintf(usage, "%s-o GAP_OPEN\n", usage);
    
    int opt = 0;
    char score_matrix_file[BUFSIZ] = "";
    char dclst_file[BUFSIZ] = "";
    char alignment_file[BUFSIZ] = "";
    int between = 0;
    int specified_line_i = -1;
    int score_for_lines = 0;
    int search_from = -1;
    int search_to = -1;
    while ((opt=getopt(argc, argv, "m:o:e:t:a:bi:lF:T:I:")) != -1) {
	switch (opt) {
	case 'm':
	    strcpy(score_matrix_file, optarg);
	    break;
	case 't':
	    strcpy(dclst_file, optarg);
	    break;
	case 'a':
	    strcpy(alignment_file, optarg);
	    break;
	case 'e':
	    GAP_EXT = atof(optarg);
	    break;
	case 'o':
	    GAP_OPEN = atof(optarg);
	    break;
	case 'b':
	    between = 1;
	    break;
	case 'i':
	    specified_line_i = atoi(optarg);
	    break;
	case 'l':
	    score_for_lines = 1;
	    break;
	case 'F':
	    search_from = atoi(optarg);
	    break;
	case 'T':
	    search_to = atoi(optarg);
	    break;
	}
    }
    argv += optind;
    argc -= optind;
    if (strcmp(score_matrix_file, "") == 0 || strcmp(dclst_file, "") == 0 || strcmp(alignment_file, "") == 0) {
	fprintf(stderr, "%s", usage);
	exit(1);
    }

    get_score_matrix(SCORE, score_matrix_file);

    /* dclst */
    int n_dclst = count_dclst(dclst_file);
    DOMAIN *dclst = (DOMAIN *)my_malloc(n_dclst * sizeof(DOMAIN));
    read_dclst(dclst_file, dclst);
    
    /* alignment */
    int n_alignment = count_alignment(alignment_file);
    char **alignment = (char **)my_malloc(n_alignment * sizeof(char *));
    char **gene = (char **)my_malloc(n_alignment * sizeof(char *));
    int n_pos = read_alignment(alignment_file, alignment, gene, n_alignment);
    fprintf(stderr, "n_seq=%d n_pos=%d ", n_alignment, n_pos);
    
    /* memory */
    int n_clusters = argc ? argc : 1;
    char **cluster_flgs = (char **)my_calloc(n_clusters, sizeof(char *));
    int c;
    for(c=0; c<n_clusters; c++){
	cluster_flgs[c] = (char *)my_calloc(n_alignment, sizeof(char));
    }
    char **domain = (char **)my_malloc(n_alignment * sizeof(char *)); /* should be int, to treat over 256 clusters at once */
    char **gap = (char **)my_malloc(n_alignment * sizeof(char *));
    int i;
    for(i=0; i<n_alignment; i++){
	domain[i] = (char *)my_calloc(n_pos, sizeof(char));
	gap[i] = (char *)my_calloc(n_pos, sizeof(char));
	map_gaps(alignment[i], gap[i], n_pos);
    }

    /* mapping */
    int k;
    for(k=0; k<n_dclst; k++){
	int cluster = argc ? get_cluster(&dclst[k], argv, argc) : 1;
	if (cluster) {
	    map_on_alignment(&dclst[k], alignment, gene, n_alignment, domain, cluster_flgs, cluster);
	}
    }
    
    /* score */
    if (score_for_lines) {
	calc_score_for_lines(alignment, n_alignment, domain, cluster_flgs, n_clusters);
    } else if (n_clusters == 2 && between) {
	double score_bet = score_between(n_alignment, alignment, domain, gap, n_clusters, cluster_flgs);
	printf("%.1f\n", score_bet);
    } else if (n_clusters == 2 && search_from >=0 && search_to >= 0 && specified_line_i >= 0) {
	if (cluster_flgs[0][specified_line_i] && cluster_flgs[1][specified_line_i]) {
	    double max_score;
	    /* int max_j = search_boundary_i(alignment, n_alignment, domain, cluster_flgs, n_clusters, specified_line_i, search_from, search_to); */
	    int max_j = search_boundary_i_local(alignment, n_alignment, domain, cluster_flgs, n_clusters, specified_line_i, search_from, search_to, &max_score);
	    /* int step = 10; */
	    /* int max_j = search_boundary_i_local_step(alignment, n_alignment, domain, cluster_flgs, n_clusters, specified_line_i, search_from, search_to, step, &max_score); */
	    /* set_boundary_i(alignment, domain, specified_line_i, max_j - step, search_from, search_to); */
	    /* max_j = search_boundary_i_local_fine(alignment, n_alignment, domain, cluster_flgs, n_clusters, specified_line_i, max(search_from, max_j - step), min(search_to, max_j + step), &max_score); */
	    printf("%d\t%.1f\n", max_j, max_score);
	} else {
	    fprintf(stderr, "cannot found boundary for gene %d\n", specified_line_i);
	}
    } else if (n_clusters == 2 && search_from >=0 && search_to >= 0) {
	double max_score;
	/* int max_pos = search_boundary(alignment, n_alignment, domain, cluster_flgs, n_clusters, search_from, search_to, &max_score); */
	int max_pos = search_boundary_local(alignment, n_alignment, domain, cluster_flgs, n_clusters, search_from, search_to, &max_score);
	printf("%d\t%.1f\n", max_pos, max_score);
    } else if (specified_line_i >= 0) {
	double score = calc_score_for_line(alignment, n_alignment, domain, cluster_flgs, n_clusters, specified_line_i);
	int n_aa = count_aa_i(alignment, domain, specified_line_i);
	printf("%f\t1\t%d\n",score, n_aa);
    } else {
	double score = calc_score(alignment, n_alignment, domain, cluster_flgs, n_clusters);
	int n_aa = count_aa(alignment, domain, n_alignment);
	printf("%.1f\t%d\t%d\n", score, n_alignment, n_aa);
    }

    /* free */

    /* free(dclst); */
    /* for(i=0; i<n_alignment; i++){ */
    /* 	free(gene[i]); */
    /* 	free(alignment[i]); */
    /* 	free(gap[i]); */
    /* 	free(domain[i]); */
    /* } */
    /* free(gene); */
    /* free(alignment); */
    /* free(gap); */
    /* free(domain); */
    /* for(c=0; c<n_clusters; c++){ */
    /* 	free(cluster_flgs[c]); */
    /* } */
    /* free(cluster_flgs); */
    
    return 0;
}

/*-----------------*/
/* Score           */
/*-----------------*/

int search_boundary_i(char **alignment, int n_alignment, char **domain, char **cluster_flgs, int n_clusters, int i, int search_from, int search_to)
{
    double max_score;
    int max_j;
    int j;
    for(j=search_from; j<=search_to; j++){
	set_boundary_i(alignment, domain, i, j, search_from, search_to);
	double score = calc_score(alignment, n_alignment, domain, cluster_flgs, n_clusters);
	if (j == search_from || score > max_score) {
	    max_score = score;
	    max_j = j;	    
	}
	/* fprintf(stderr, "j=%d score=%.1f\n", j, score); */
    }
    return max_j;
}

int search_boundary_i_local(char **alignment, int n_alignment, char **domain, char **cluster_flgs, int n_clusters, int i, int search_from, int search_to, double *max_score)
{
    int j = search_from;
    set_boundary_i(alignment, domain, i, j, search_from, search_to);
    double score = calc_score(alignment, n_alignment, domain, cluster_flgs, n_clusters);
        
    int max_j = j;
    *max_score = score;
    for(j=search_from+1; j<=search_to; j++){
	double score_before = sum_of_pairs_for_line_local(alignment, domain, cluster_flgs, n_alignment, n_clusters, i, j-1, j);
	int i2;
	int n_gaps_before = 0;
	n_gaps_before += count_gaps_local(alignment, domain, cluster_flgs, i, n_alignment, max(0,j-2), j);
	for(i2=0; i2<n_alignment; i2++){
	    if (i != i2) {
		n_gaps_before += count_gaps_local_pair(alignment, domain, cluster_flgs, i2, i, n_alignment, max(0,j-2), j);
	    }
	}
	set_boundary_i(alignment, domain, i, j, j-1, j);
	double score_after = sum_of_pairs_for_line_local(alignment, domain, cluster_flgs, n_alignment, n_clusters, i, j-1, j);
	int n_gaps_after = 0;
	n_gaps_after += count_gaps_local(alignment, domain, cluster_flgs, i, n_alignment, max(0,j-2), j);
	for(i2=0; i2<n_alignment; i2++){
	    if (i != i2) {
		n_gaps_after += count_gaps_local_pair(alignment, domain, cluster_flgs, i2, i, n_alignment, max(0,j-2), j);
	    }
	}
	score += (score_after - score_before) - (n_gaps_after - n_gaps_before) * GAP_OPEN;
	if (score > *max_score) {
	    *max_score = score;
	    max_j = j;	    
	}
    }
    return max_j;
}

int search_boundary_i_local_fine(char **alignment, int n_alignment, char **domain, char **cluster_flgs, int n_clusters, int i, int search_from, int search_to, double *max_score)
{
    int j = search_from;
    double score = calc_score(alignment, n_alignment, domain, cluster_flgs, n_clusters);
        
    int max_j = j;
    *max_score = score;
    for(j=search_from+1; j<=search_to; j++){
	double score_before = sum_of_pairs_for_line_local(alignment, domain, cluster_flgs, n_alignment, n_clusters, i, j-1, j);
	int i2;
	int n_gaps_before = 0;
	n_gaps_before += count_gaps_local(alignment, domain, cluster_flgs, i, n_alignment, max(0,j-2), j);
	for(i2=0; i2<n_alignment; i2++){
	    if (i != i2) {
		n_gaps_before += count_gaps_local_pair(alignment, domain, cluster_flgs, i2, i, n_alignment, max(0,j-2), j);
	    }
	}
	set_boundary_i(alignment, domain, i, j, j-1, j);
	double score_after = sum_of_pairs_for_line_local(alignment, domain, cluster_flgs, n_alignment, n_clusters, i, j-1, j);
	int n_gaps_after = 0;
	n_gaps_after += count_gaps_local(alignment, domain, cluster_flgs, i, n_alignment, max(0,j-2), j);
	for(i2=0; i2<n_alignment; i2++){
	    if (i != i2) {
		n_gaps_after += count_gaps_local_pair(alignment, domain, cluster_flgs, i2, i, n_alignment, max(0,j-2), j);
	    }
	}
	score += (score_after - score_before) - (n_gaps_after - n_gaps_before) * GAP_OPEN;
	if (score > *max_score) {
	    *max_score = score;
	    max_j = j;	    
	}
    }
    return max_j;
}

int search_boundary_i_local_step(char **alignment, int n_alignment, char **domain, char **cluster_flgs, int n_clusters, int i, int search_from, int search_to, int step, double *max_score)
{
    int j = search_from;
    set_boundary_i(alignment, domain, i, j, search_from, search_to);
    double score = calc_score(alignment, n_alignment, domain, cluster_flgs, n_clusters);
        
    int max_j = j;
    *max_score = score;

    j+=step;
    while(j<=search_to){
	double score_before = sum_of_pairs_for_line_local(alignment, domain, cluster_flgs, n_alignment, n_clusters, i, j-step, j);
	int i2;
	int n_gaps_before = 0;
	n_gaps_before += count_gaps_local(alignment, domain, cluster_flgs, i, n_alignment, max(0,j-step-1), j);
	for(i2=0; i2<n_alignment; i2++){
	    if (i != i2) {
		n_gaps_before += count_gaps_local_pair(alignment, domain, cluster_flgs, i2, i, n_alignment, max(0,j-step-1), j);
	    }
	}
	set_boundary_i(alignment, domain, i, j, j-step, j);
	double score_after = sum_of_pairs_for_line_local(alignment, domain, cluster_flgs, n_alignment, n_clusters, i, j-step, j);
	int n_gaps_after = 0;
	n_gaps_after += count_gaps_local(alignment, domain, cluster_flgs, i, n_alignment, max(0,j-step-1), j);
	for(i2=0; i2<n_alignment; i2++){
	    if (i != i2) {
		n_gaps_after += count_gaps_local_pair(alignment, domain, cluster_flgs, i2, i, n_alignment, max(0,j-step-1), j);
	    }
	}
	score += (score_after - score_before) - (n_gaps_after - n_gaps_before) * GAP_OPEN;
	if (score > *max_score) {
	    *max_score = score;
	    max_j = j;	    
	}
	j+=step;
    }
    return max_j;
}

int search_boundary_local (char **alignment, int n_alignment, char **domain, char **cluster_flgs, int n_clusters, int search_from, int search_to, double *max_score)
{

    int j = search_from;
    set_boundary(alignment, n_alignment, domain, cluster_flgs, search_from, search_to, j);
    double score = calc_score(alignment, n_alignment, domain, cluster_flgs, n_clusters);
    fprintf(stderr, "j=%d score=%.1f\n", j, score);

    int max_j = j;
    *max_score = score;

    int n_pos = strlen(alignment[0]);
    for(j=search_from+1; j<=search_to; j++){
	int i;
	int n_gaps_before = 0;
	for(i=0; i<n_alignment; i++){
	    n_gaps_before += count_gaps_local_with_flg(alignment, domain, cluster_flgs, i, n_alignment, max(0,j-2), j);
	}
	for(i=0; i<n_alignment; i++){
	    if (cluster_flgs[0][i] && cluster_flgs[1][i]) {
		double score_before = sum_of_pairs_for_line_local(alignment, domain, cluster_flgs, n_alignment, n_clusters, i, j-1, j);
		set_boundary_i(alignment, domain, i, j, j-1, j);
		double score_after = sum_of_pairs_for_line_local(alignment, domain, cluster_flgs, n_alignment, n_clusters, i, j-1, j);
		score += (score_after - score_before);
	    }
	}
	int n_gaps_after = 0;
	for(i=0; i<n_alignment; i++){
	    n_gaps_after += count_gaps_local_with_flg(alignment, domain, cluster_flgs, i, n_alignment, max(0,j-2), j);
	}
	score += - (n_gaps_after - n_gaps_before) * GAP_OPEN;
	if (score > *max_score) {
	    *max_score = score;
	    max_j = j;	    
	}
	fprintf(stderr, "j=%d score=%.1f\n", j, score);
    }
    return max_j;
}

int search_boundary (char **alignment, int n_alignment, char **domain, char **cluster_flgs, int n_clusters, int search_from, int search_to, double *max_score)
{
    int max_j;
    int j;
    for(j=search_from; j<=search_to; j++){
	set_boundary(alignment, n_alignment, domain, cluster_flgs, search_from, search_to, j);
	double score = calc_score(alignment, n_alignment, domain, cluster_flgs, n_clusters);
	if (j == search_from || score > *max_score) {
	    *max_score = score;
	    max_j = j;	    
	}
	fprintf(stderr, "j=%d score=%.1f\n", j, score);
    }
    return max_j;
}

double calc_score_for_line (char **alignment, int n_alignment, char **domain, char **cluster_flgs, int n_clusters, int i)
{
    int n_pos = strlen(alignment[0]);
    double sum = sum_of_pairs_for_line_local(alignment, domain, cluster_flgs, n_alignment, n_clusters, i, 0, n_pos-1);
    int n_gaps = 0;
    int cluster;
    for(cluster=1; cluster<=n_clusters; cluster++){

	n_gaps += count_gaps(alignment, domain, cluster_flgs, i, n_alignment, cluster);

	int i2;
	for(i2=0; i2<n_alignment; i2++){
	    if (i != i2) {
		n_gaps += count_gaps_pair(alignment, domain, cluster_flgs, i2, i, n_alignment, cluster);
	    }
	}

    }
    return sum/2 - n_gaps * GAP_OPEN;
}

double calc_score_for_lines (char **alignment, int n_alignment, char **domain, char **cluster_flgs, int n_clusters)
{
    int n_pos = strlen(alignment[0]);
    int i;
    for(i=0; i<n_alignment; i++){
	double sum = sum_of_pairs_for_line_local(alignment, domain, cluster_flgs, n_alignment, n_clusters, i, 0, n_pos-1);
	double score_for_line = sum/2;
	int n_aa = count_aa_i(alignment, domain, i);
	if (i != 0) {
	    printf(" ");
	}
	printf("%f",score_for_line/n_alignment/n_aa);
    }
    printf("\n");
}

double calc_score(char **alignment, int n_alignment, char **domain, char **cluster_flgs, int n_clusters) 
{
    double score = 0;
    int cluster;
    double sum;
    int n_gaps;
    int i, i2, j;
    int n_pos = strlen(alignment[0]);
    for(cluster=1; cluster<=n_clusters; cluster++){
	sum = 0;
	n_gaps = 0;
	for(i=0; i<n_alignment; i++){
	    if (cluster_flgs[cluster-1][i]) {
		for(i2=i+1; i2<n_alignment; i2++){
		    if (cluster_flgs[cluster-1][i2]) {
			for(j=0; j<n_pos; j++){
			    if ((alignment[i][j] != '-' && domain[i][j] == cluster) && (alignment[i2][j] != '-' && domain[i2][j] == cluster)) {
				sum += SCORE[alignment[i][j]][alignment[i2][j]];
			    } else if ((alignment[i][j] != '-' && domain[i][j] == cluster) || (alignment[i2][j] != '-' && domain[i2][j] == cluster)) {
				sum -= GAP_EXT;
			    }
			}
		    }
		}
	    }
	    n_gaps += count_gaps(alignment, domain, cluster_flgs, i, n_alignment, cluster);
	}
	/* fprintf(stderr, "[%d] %.1f - %d * %.1f\n", cluster, sum, n_gaps, GAP_OPEN); */
	score += sum - n_gaps * GAP_OPEN;
    }
    return score;
}

double score_between (int n_alignment, char **alignment, char **domain, char **gap, int n_clusters, char **cluster_flgs)
{
    double score_max = 0;
    int max_flg = 0;
    int i, i2;
    for(i=0; i<n_alignment; i++){
	for(i2=0; i2<n_alignment; i2++){
	    if (i != i2 && cluster_flgs[0][i] && cluster_flgs[1][i2]) {
		double score = sum_of_pairs_between(alignment, domain, gap, cluster_flgs, n_alignment, n_clusters, i, i2);
		int n_open = 0;
		n_open += count_gaps_between(alignment, domain, cluster_flgs, i, i2, n_alignment, n_clusters, 1, 2);
		n_open += count_gaps_between(alignment, domain, cluster_flgs, i2, i, n_alignment, n_clusters, 2, 1);
		score -= n_open * GAP_OPEN;
		if (! max_flg || score > score_max) {
		    score_max = score;
		    max_flg = 1;
		}
	    }
	}
    }
    return score_max;
}

double sum_of_pairs_for_line_local(char **a, char **d, char **cluster_flgs, int n, int n_clusters, int i, int from, int to)
{
    double sum = 0;
    int i2, j, cluster;
    for(cluster=1; cluster<=n_clusters; cluster++){
	if (cluster_flgs[cluster-1][i]) {
	    for(i2=0; i2<n; i2++){
		if (i != i2) {
		    if (cluster_flgs[cluster-1][i2]) {
			for(j=from; j<=to; j++){
			    if ((a[i][j] != '-' && d[i][j] == cluster) && (a[i2][j] != '-' && d[i2][j] == cluster)) {
				sum += SCORE[a[i][j]][a[i2][j]];
			    } else if ((a[i][j] != '-' && d[i][j] == cluster) || (a[i2][j] != '-' && d[i2][j] == cluster)) {
				sum -= GAP_EXT;
			    }
			}
		    }
		}
	    }
	}
    }

    return sum;
}

double sum_of_pairs_between(char **a, char **d, char **g, char **cluster_flgs, int n, int n_clusters, int i, int i2)
{
    double sum = 0;
    int n_pos = strlen(a[0]);
    int j;
    for(j=0; j<n_pos; j++){
	if ((a[i][j] != '-' && d[i][j] == 1) && (a[i2][j] != '-' && d[i2][j] == 2)) {
	    sum += SCORE[a[i][j]][a[i2][j]];
	} else if ((a[i][j] != '-' && d[i][j] == 1)                     && (g[i2][j] == 1   && d[i2][j] != 2 || a[i2][j] != '-')) {
	} else if ((g[i][j] == 1   && d[i][j] != 1  ||  a[i][j] != '-') && (a[i2][j] != '-' && d[i2][j] == 2)                   ) {
	} else if ((a[i][j] != '-' && d[i][j] == 1) || (a[i2][j] != '-' && d[i2][j] == 2)) {
	    sum -= GAP_EXT;
	}
    }

    return sum;
}

int count_gaps(char **alignment, char **domain, char **cluster_flgs, int i, int n, int cluster)
{
    int n_open = 0;

    int m = strlen(alignment[0]);
    if (cluster_flgs[cluster-1][i]) {
	int start;
	int end;
	int j = 0;
	while (j<m) {
	    if (alignment[i][j] == '-' || domain[i][j] != cluster) {
		start = j;
		while (j<m && (alignment[i][j] == '-' || domain[i][j] != cluster)) {
		    j++;
		}
		end = j - 1;
		n_open += eval_gap(alignment, domain, cluster_flgs, i, n, start, end, cluster);
	    }
	    j++;
	}
    }    
    
    return n_open;
}

int count_gaps_pair(char **alignment, char **domain, char **cluster_flgs, int i, int i2, int n, int cluster)
{
    int n_open = 0;

    int m = strlen(alignment[0]);
    if (cluster_flgs[cluster-1][i] && cluster_flgs[cluster-1][i2]) {
	int start;
	int end;
	int j = 0;
	while (j<m) {
	    if (alignment[i][j] == '-' || domain[i][j] != cluster) {
		start = j;
		while (j<m && (alignment[i][j] == '-' || domain[i][j] != cluster)) {
		    j++;
		}
		end = j - 1;
		n_open += find_aa(alignment[i2], domain[i2], start, end, cluster);
	    }
	    j++;
	}
    }
    
    return n_open;
}

int count_gaps_local(char **alignment, char **domain, char **cluster_flgs, int i, int n, int from, int to)
{
    int n_open = 0;

    int m = strlen(alignment[0]);
    int cluster;
    for(cluster=1; cluster<=2; cluster++){
	if (cluster_flgs[cluster-1][i]) {
	    int start;
	    int end;
	    int j = from;
	    if (alignment[i][j] == '-' || domain[i][j] != cluster) {
		while (j>=0 && (alignment[i][j] == '-' || domain[i][j] != cluster)) {
		    j--;
		}
		start = j + 1;
		j = from;
		while (j<m && (alignment[i][j] == '-' || domain[i][j] != cluster)) {
		    j++;
		}
		end = j - 1;
		n_open += eval_gap(alignment, domain, cluster_flgs, i, n, start, end, cluster);
	    }
	    while (j<=to) {
		if (alignment[i][j] == '-' || domain[i][j] != cluster) {
		    start = j;
		    while (j<m && (alignment[i][j] == '-' || domain[i][j] != cluster)) {
			j++;
		    }
		    end = j - 1;
		    n_open += eval_gap(alignment, domain, cluster_flgs, i, n, start, end, cluster);
		}
		j++;
	    }
	}    
    }

    return n_open;
}

int count_gaps_local_pair(char **alignment, char **domain, char **cluster_flgs, int i, int i2, int n, int from, int to)
{
    int n_open = 0;

    int m = strlen(alignment[0]);
    int cluster;
    for(cluster=1; cluster<=2; cluster++){
	if (cluster_flgs[cluster-1][i] && cluster_flgs[cluster-1][i2]) {
	    int start;
	    int end;
	    int j = from;
	    if (alignment[i][j] == '-' || domain[i][j] != cluster) {
		while (j>=0 && (alignment[i][j] == '-' || domain[i][j] != cluster)) {
		    j--;
		}
		start = j + 1;
		j = from;
		while (j<m && (alignment[i][j] == '-' || domain[i][j] != cluster)) {
		    j++;
		}
		end = j - 1;
		n_open += find_aa(alignment[i2], domain[i2], start, end, cluster);
	    }
	    while (j<=to) {
		if (alignment[i][j] == '-' || domain[i][j] != cluster) {
		    start = j;
		    while (j<m && (alignment[i][j] == '-' || domain[i][j] != cluster)) {
			j++;
		    }
		    end = j - 1;
		    n_open += find_aa(alignment[i2], domain[i2], start, end, cluster);
		}
		j++;
	    }
	}    
    }

    return n_open;
}

int count_gaps_local_with_flg(char **alignment, char **domain, char **cluster_flgs, int i, int n, int from, int to)
{
    int n_open = 0;

    int m = strlen(alignment[0]);
    int cluster;
    for(cluster=1; cluster<=2; cluster++){
	if (cluster_flgs[cluster-1][i]) {
	    int start;
	    int end;
	    int j = from;
	    if (alignment[i][j] == '-' || domain[i][j] != cluster) {
		while (j>=0 && (alignment[i][j] == '-' || domain[i][j] != cluster)) {
		    j--;
		}
		start = j + 1;
		j = from;
		while (j<m && (alignment[i][j] == '-' || domain[i][j] != cluster)) {
		    j++;
		}
		end = j - 1;
		n_open += eval_gap_with_flg(alignment, domain, cluster_flgs, i, n, start, end, cluster);
	    }
	    while (j<=to) {
		if (alignment[i][j] == '-' || domain[i][j] != cluster) {
		    start = j;
		    while (j<m && (alignment[i][j] == '-' || domain[i][j] != cluster)) {
			j++;
		    }
		    end = j - 1;
		    n_open += eval_gap_with_flg(alignment, domain, cluster_flgs, i, n, start, end, cluster);
		}
		j++;
	    }
	}    
    }

    return n_open;
}

int count_gaps_between(char **alignment, char **domain, char **cluster_flgs, int i, int i2, int n, int n_clusters, int cluster1, int cluster2)
{
    /* does not use cluster_flgs ? */
    int n_open = 0;

    int m = strlen(alignment[0]);
    int start;
    int end;
    int j = 0;
    while (j<m) {
	if (alignment[i][j] == '-' || domain[i][j] != cluster1) {
	    start = j;
	    while (j<m && (alignment[i][j] == '-' || domain[i][j] != cluster1)) {
		j++;
	    }
	    end = j - 1;
	    n_open += find_aa(alignment[i2], domain[i2], start, end, cluster2);
	}
	j++;
    }
    
    return n_open;
}

int eval_gap(char **alignment, char **domain, char **cluster_flgs, int i, int n, int start, int end, int cluster)
{
    int n_open = 0;
    int i2;
    for(i2=0; i2<n; i2++){
	if (i != i2 && cluster_flgs[cluster-1][i2]) {
	    n_open += find_aa(alignment[i2], domain[i2], start, end, cluster);
	}
    }
    return n_open;
}

int eval_gap_with_flg(char **alignment, char **domain, char **cluster_flgs, int i, int n, int start, int end, int cluster)
{
    int n_open = 0;
    int i2;
    for(i2=0; i2<n; i2++){
	if (cluster_flgs[0][i] && cluster_flgs[1][i] ||
	    cluster_flgs[0][i2] && cluster_flgs[1][i2]
	    ) {
	    if (i != i2 && cluster_flgs[cluster-1][i2]) {
		n_open += find_aa(alignment[i2], domain[i2], start, end, cluster);
	    }
	}
    }
    return n_open;
}

int find_aa(char alignment[], char domain[], int start, int end, int cluster)
{
    int j;
    for(j=start; j<=end; j++){
	if (alignment[j] != '-' && domain[j] == cluster) {
	    return 1;
	}
    }
    return 0;
}

/*-----------------*/
/* Functions       */
/*-----------------*/

void set_boundary(char **alignment, int n_alignment, char **domain, char **cluster_flgs, int search_from, int search_to, int j_boundary)
{
    int i;
    for(i=0; i<n_alignment; i++){
	if (cluster_flgs[0][i] && cluster_flgs[1][i]) {
	    set_boundary_i(alignment, domain, i, j_boundary, search_from, search_to);
	}
    }
}

/* bug ? */
void set_boundary_i(char **alignment, char **domain, int i, int j_boundary, int search_from, int search_to)
{
    int j;
    j = j_boundary - 1;
    while (alignment[i][j] == '-') {
	domain[i][j] = 0;
	j --;
    }
    while (j >= search_from) {
	domain[i][j] = 1;
	j--;
    }
    j = j_boundary;
    while (alignment[i][j] == '-') {
	domain[i][j] = 0;
	j ++;
    }
    while (j <= search_to) {
	domain[i][j] = 2;
	j++;   
    }
}

void set_boundary_i_slow(char **alignment, char **domain, int i, int j_boundary, int search_from, int search_to)
{
    int j;

    j = j_boundary - 1;
    while (alignment[i][j] == '-') {
	domain[i][j] = 0;
	j --;
    }
    while (j >= search_from) {
	domain[i][j] = 1;
	j--;
    }
    while (alignment[i][j] == '-') {
	domain[i][j] = 1;
	j--;
    }

    j = j_boundary;
    while (alignment[i][j] == '-') {
	domain[i][j] = 0;
	j ++;
    }
    while (j <= search_to) {
	domain[i][j] = 2;
	j++;   
    }
    while (alignment[i][j] == '-') {
	domain[i][j] = 2;
	j++;
    }
}

void map_on_alignment(DOMAIN *a_domain, char **alignment, char **gene, int n_alignment, char **domain, char **cluster_flgs, int cluster)
{
    int n_pos = strlen(alignment[0]);
    int i,j;
    int n_mapped = 0;
    for(i=0; i<n_alignment; i++){
	if (strcmp(a_domain->gene, gene[i]) == 0) {
	    int pos = 0;
	    for(j=0; j<n_pos; j++){
		if (alignment[i][j] != '-') {
		    pos++;
		}
		if (pos >= a_domain->begin && pos < a_domain->end ||
		    alignment[i][j] != '-' && pos == a_domain->end
		    ) {
		    domain[i][j] = cluster;
		}
	    }
	    n_mapped ++;
	    cluster_flgs[cluster-1][i] = 1;
	}
    }
    if (n_mapped != 1) {
	fprintf(stderr, "WARNING: domain info of %s mapped to %d sequences in the alignment\n", a_domain->gene, n_mapped);
	/* exit(1); */
    }
}

void map_gaps (char *alignment, char *gap, int n_pos)
{
    /* 1: internal gap, 2: terminal gap */
    int j = 0;
    while (alignment[j] == '-') {
	gap[j] = 2;
	j++;	
    }
    int j1 = j;
    j = n_pos - 1;
    while (alignment[j] == '-') {
	gap[j] = 2;
	j--;
    }
    int j2 = j;
    for(j=j1; j<=j2; j++){
	if (alignment[j] == '-') {
	    gap[j] = 1;
	}
    }
}

int count_aa (char **a, char **d, int n)
{
    int n_aa = 0;
    int i;
    for(i=0; i<n; i++){
	n_aa += count_aa_i(a, d, i);	
    }
    return n_aa;
}

int count_aa_i (char **a, char **d, int i)
{
    int n_aa = 0;
    int m = strlen(a[0]);
    int j;
    for(j=0; j<m; j++){
	if (a[i][j] != '-') {
	    if (d[i][j] != 0) {
		n_aa ++;
	    }
	}
    }
    return n_aa;
}

int get_cluster (DOMAIN *a_domain, char **clusters, int n_clusters)
{
    int cluster = 0;

    int c;
    for(c=0; c<n_clusters; c++){
	if (strcmp(a_domain->cluster, clusters[c]) == 0) {
	    if (cluster) {
		fprintf(stderr, "cluster %s duplicated\n", a_domain->cluster);
		exit(10);
	    }
	    cluster = c + 1;
	}
    }

    return cluster;
}

/*-----------------*/
/* Read Functions  */
/*-----------------*/

/* effective enough ? error handling is enough ? */
void get_score_matrix(int score[][256], char *file)
{
    char key[256];
    int mat[256][256];
    int n = read_matrix(file, key, mat);
    
    int i, j;
    for(i=0; i<256; i++){
	for(j=0; j<256; j++){
	    score[i][j] = 0;
	}
    }
    for(i=0; i<n; i++){
	for(j=0; j<n; j++){
	    score[key[i]][key[j]] = mat[i][j];
	}
    }
    
}

int read_matrix (char *file, char key[], int mat[][256])
{
    FILE *fp;
    if ((fp = fopen(file, "r")) == NULL) {
	fprintf(stderr, "can't open %s\n", file);
	exit(2);
    }
    
    char buf[BUFSIZ] = "";
    int array[BUFSIZ];
    
    int n;
    int i = 0;
    while (fgets(buf, BUFSIZ, fp)) {
	chomp(buf);
	if (buf[0] == '#') {
	} else if (buf[0] == ' ') {
	    n = split_char(&buf[2], key);
	} else if (buf[0] != ' ' && buf[1] == ' ') {
	    n = split_int(&buf[2], mat[i]);
	    i++;
	} else {
	    exit(3);
	}
    }
    
    return i;
}

int read_alignment (char *file, char **alignment, char **gene, int n_alignment)
{
    FILE *fp;
    if ((fp = fopen(file, "r")) == NULL) {
	fprintf(stderr, "can't open %s\n", file);
	exit(2);
    }

    char buf[BUFSIZ] = "";
    int bufsiz = BUFSIZ;
    int i = -1;
    while (fgets(buf, BUFSIZ, fp)) {
	chomp(buf);
	if (buf[0] == '>') {
	    i++;
	    if (i >= n_alignment) {
		fprintf(stderr, "exceeded n_alignment = %d\n", n_alignment);
		exit(1);
	    }
	    alignment[i] = (char *)my_malloc(bufsiz * sizeof(char));
	    alignment[i][0] = '\0';
	    gene[i] = (char *)my_malloc(strlen(buf) * sizeof(char));
	    strcpy(gene[i], &buf[1]); /* remove the first char '>' */
	} else {
	    if ((int)strlen(alignment[i]) + (int)strlen(buf) >= bufsiz) {
		/* fprintf(stderr, "alignment_len: %d + %d >= %d, ", (int)strlen(alignment[i]), (int)strlen(buf), bufsiz); */
		bufsiz *= 2;
		/* fprintf(stderr, "realloc %d\n", bufsiz); */
		alignment[i] = (char *)my_realloc(alignment[i], bufsiz * sizeof(char));
	    }
	    strcat(alignment[i], buf);
	}
    }

    fclose(fp);
    
    int n_pos = strlen(alignment[0]);
    return n_pos;
}

int count_alignment (char *file)
{
    FILE *fp;
    if ((fp = fopen(file, "r")) == NULL) {
	fprintf(stderr, "can't open %s\n", file);
	exit(2);
    }
    
    char buf[BUFSIZ] = "";
    int i = 0;
    while (fgets(buf, BUFSIZ, fp)) {
	if (buf[0] == '>') {
	    i++;
	}
    }
    fclose(fp);

    return i;
}

int count_dclst (char *file)
{
    FILE *fp;
    if ((fp = fopen(file, "r")) == NULL) {
	fprintf(stderr, "can't open %s\n", file);
	exit(2);
    }
    
    char buf[BUFSIZ] = "";
    int i = 0;
    while (fgets(buf, BUFSIZ, fp)) {
	i++;
    }

    return i;
}

int read_dclst (char *file, DOMAIN dclst[])
{
    FILE *fp;
    if ((fp = fopen(file, "r")) == NULL) {
	fprintf(stderr, "can't open %s\n", file);
	exit(2);
    }
    
    char buf[BUFSIZ] = "";
    int i = 0;
    while (fgets(buf, BUFSIZ, fp)) {
	chomp(buf);
	split_line_to_domain(buf, &dclst[i]);
	i++;
    }
    
    return i;
}

/*-------------------*/
/* General Functions */
/*-------------------*/

void *my_calloc (size_t n, size_t size)
{
    void *ret = calloc(n, size);
    if (ret == NULL) {
	fprintf(stderr, "can't allocate memory\n");
	exit(EXIT_FAILURE);
    }
    return ret;
}

void *my_realloc (void *ptr, size_t size)
{
    void *ret = realloc(ptr, size);
    if (ret == NULL) {
	fprintf(stderr, "can't allocate memory\n");
	exit(EXIT_FAILURE);
    }
    return ret;
}

void *my_malloc (size_t size)
{
    void *ret = malloc(size);
    if (ret == NULL) {
	fprintf(stderr, "can't allocate memory\n");
	exit(EXIT_FAILURE);
    }
    return ret;
}

void chomp(char *s)
{
    int n = strlen(s);
    if (s[n-1] == '\n') {
	s[n-1] = '\0';
    }
}

int split_int(char *line, int array[])
{
    char *str = strtok(line, " ");
    int i = 0;
    while (str != NULL) {
        array[i++] = atoi(str);
        str = strtok(NULL, " ");
    }
    return i;
}

int split_csv_to_int(char *line, int *array)
{
    char *str = strtok(line, ",");
    int i = 0;
    while (str != NULL) {
        array[i++] = atoi(str);
        str = strtok(NULL, ",");
    }
    return i;
}

int split_char(char *line, char array[])
{
    char *str = strtok(line, " ");
    int i = 0;
    while (str != NULL) {
        array[i++] = str[0];
        str = strtok(NULL, " ");
    }
    return i;
}

int split_line (char *line, char f[][BUFSIZ])
{
    char *str;
    int i = 0;
    str = strtok(line, " \t");
    while (str != NULL) {
	if (strlen(str) >= BUFSIZ) {
	    fprintf(stderr, "%s\n", line);
	    exit(3);
	}
	strcpy(f[i], str);
	i++;
	str = strtok(NULL, " \t");
    }
    return i;
}

int split_line_to_domain(char *line, DOMAIN *dclst)
{
    char f[5][BUFSIZ];
    int n = split_line(line, f);
    if (n == 5) {
	strcpy(dclst->cluster, f[0]);
	strcpy(dclst->gene, f[1]);
	dclst->domain = atoi(f[2]);
	dclst->begin = atoi(f[3]);
	dclst->end = atoi(f[4]);
    } else if (n == 4) {
	strcpy(dclst->cluster, f[0]);
	strcpy(dclst->gene, f[1]);
	dclst->domain = 0;
	dclst->begin = atoi(f[2]);
	dclst->end = atoi(f[3]);
    } else {
	fprintf(stderr, "%s\n", line);
	exit(3);
    }
}

/* int split_line_to_domain(char *line, DOMAIN *dclst) */
/* { */
/*     char *str = strtok(line, " \t"); */
/*     strcpy(dclst->cluster, str); */

/*     str = strtok(NULL, " \t"); */
/*     strcpy(dclst->gene, str); */

/*     str = strtok(NULL, " \t"); */
/*     dclst->domain = atoi(str); */

/*     str = strtok(NULL, " \t"); */
/*     dclst->begin = atoi(str); */

/*     str = strtok(NULL, " \t"); */
/*     dclst->end = atoi(str); */

/*     str = strtok(NULL, " \t"); */
/*     if (str != NULL) { */
/* 	fprintf(stderr, "%s\n", line); */
/* 	exit(3); */
/*     } */
/* } */

/*------------------*/
/* Print Functions  */
/*------------------*/

void print_char(char array[], int n)
{
    int i;
    /* printf("%c", array[0]); */
    printf("%d", array[0]);
    for(i=1; i<n; i++){
	/* printf("\t%c ", array[i]); */
	printf("%d", array[i]);
    }
    printf("\n");
}

void print_int(int array[], int n)
{
    int i;
    printf("%d", array[0]);
    for(i=1; i<n; i++){
	printf("\t%d ", array[i]);
    }
    printf("\n");
}

void print_mat(int mat[][256], int n)
{
    int i;
    for(i=0; i<n; i++){
	print_int(mat[i], n);
    }
}

int max (int x, int y)
{
    return (x > y) ? x : y;
}

int min (int x, int y)
{
    return (x < y) ? x : y;
}
