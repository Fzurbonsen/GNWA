/*

    gnwa.h
    This file contains the definitions for the Graph-Needleman-Wunsch-Algorithm structs and functions.
    The design of this library is heavily inspired by gssw.
    Author: Frederic zur Bonsen <fzurbonsen@student.ethz.ch>

*/

#ifndef GNWA_H
#define GNWA_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>


typedef struct {
    char type;
    uint32_t length;
} gnwa_cigar_element_t;


typedef struct {
    int32_t length;
    gnwa_cigar_element_t* elements;
} gnwa_cigar_t;


typedef struct _gnwa_node_t gnwa_node_t;
typedef struct _gnwa_node_t {
    uint64_t id; // id of the node
    const char* seq; // sequence
    int8_t* num; // numerical conversion of seq
    int32_t len; // length of the sequence
    gnwa_node_t** prev; // previous nodes
    int32_t count_prev; // number of previous nodes
    gnwa_node_t** next; // next nodes
    int32_t count_next; // number of next nodes
} _gnwa_node_t;


typedef struct {
    uint32_t size; // size of the graph
    gnwa_node_t* max_node; // node at the begining of the nodes
    gnwa_node_t** nodes; // nodes contained in the graph
} gnwa_graph_t;


typedef struct {
    uint32_t len; // length of the path
    gnwa_node_t** nodes; // nodes in the path in order
} gnwa_path_t;


typedef struct {
    int32_t n_paths;
    int32_t counter;
    gnwa_path_t** paths;
} gnwa_graph_paths_t;


typedef struct {
    const char* read; // read
    gnwa_graph_t* graph; // graph
    int32_t score; // score of the alignment
    gnwa_cigar_t cigar; // CIGAR corresponding to the alignment
    gnwa_path_t* path; // Path of the alignment
} gnwa_alignment_t;



#ifdef __cplusplus
extern "C" {
#endif // __cplusplus


int32_t gnwa_max(int32_t a, int32_t b, int32_t c);

int8_t* gnwa_create_num(const char* seq,
                        const int32_t len,
                        const int8_t* nt_table);
void gnwa_num_print(FILE* file, int8_t* num, int32_t len);

int8_t* gnwa_create_score_matrix(int32_t match, int32_t mismatch);
int8_t* gnwa_create_nt_table(void);

void gnwa_cigar_print(FILE* file, gnwa_cigar_t* cigar);
void gnwa_cigar_destroy(gnwa_cigar_t* cigar);

gnwa_node_t* gnwa_node_create(const uint64_t id,
                                const char* seq,
                                int8_t* nt_table,
                                int8_t* score_matrix);
void gnwa_node_destroy(gnwa_node_t* node);
void gnwa_node_add_prev(gnwa_node_t* node, gnwa_node_t* prev);
void gnwa_node_add_next(gnwa_node_t* node, gnwa_node_t* next);
void gnwa_node_del_prev(gnwa_node_t* node, gnwa_node_t* prev);
void gnwa_node_del_next(gnwa_node_t* node, gnwa_node_t* next);
void gnwa_node_add_edge(gnwa_node_t* start, gnwa_node_t* end);
void gnwa_node_del_edge(gnwa_node_t* start, gnwa_node_t* end);

gnwa_path_t* gnwa_path_create(gnwa_node_t** nodes, uint32_t len);
void gnwa_path_destroy(gnwa_path_t* path);
char* gnwa_path_get_sequence(gnwa_path_t* path);
int32_t gnwa_path_get_num_sequence(gnwa_path_t* path, int8_t** num);
void gnwa_path_print(FILE* file, gnwa_path_t* path);

gnwa_graph_t* gnwa_graph_create(uint32_t size);
void gnwa_graph_destroy(gnwa_graph_t* graph);
uint32_t gnwa_graph_add_node(gnwa_graph_t* graph, gnwa_node_t* node);
gnwa_node_t* gnwa_graph_find_max_node(gnwa_node_t** nodes, int32_t n_nodes);
uint32_t gnwa_graph_get_path_count(gnwa_graph_t* graph);
gnwa_graph_paths_t* gnwa_graph_get_paths(gnwa_graph_t* graph);
void gnwa_grpah_paths_destroy(gnwa_graph_paths_t* paths);
void gnwa_graph_print(FILE* file, gnwa_graph_t* graph);


gnwa_alignment_t* gnwa_alignment_create(const char* read,
                                        gnwa_graph_t* graph,
                                        int32_t score,
                                        gnwa_path_t* path);
void gnwa_alignment_destroy(gnwa_alignment_t* alignment);



gnwa_alignment_t* gnwa_align(const int8_t* num_read,
                                const int32_t read_len,
                                const int8_t* num_seq,
                                const int32_t seq_len,
                                gnwa_alignment_t* alignment,
                                int8_t* nt_table,
                                int8_t* score_matrix,
                                int8_t gap_open,
                                int8_t gap_extend);

gnwa_alignment_t* gnwa_path_align(const char* read,
                                    gnwa_graph_t* graph,
                                    gnwa_path_t* path,
                                    int8_t* nt_table,
                                    int8_t* score_matrix,
                                    uint8_t gap_open,
                                    uint8_t gap_extend);

gnwa_alignment_t* gnwa_path_build_and_align(const char* read,
                                            gnwa_node_t** nodes,
                                            gnwa_alignment_t* alignment,
                                            gnwa_graph_t* graph,
                                            int32_t len,
                                            int8_t* nt_table,
                                            int8_t* score_matrix,
                                            uint8_t gap_open,
                                            uint8_t gap_extend);

gnwa_alignment_t* gnwa_graph_align(const char* read,
                                    gnwa_graph_t* graph,
                                    int8_t* nt_table,
                                    int8_t* score_matrix,
                                    uint8_t gap_open,
                                    uint8_t gap_extend);



#ifdef __cplusplus
}
#endif	// __cplusplus

#endif // GNWA_H