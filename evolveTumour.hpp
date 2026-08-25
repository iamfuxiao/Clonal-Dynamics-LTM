/*
    File: evolveTumour.hpp
    Model: particleCell -- lineage tracing model (LTM)
    Created: 10 May, 2018 (XF)
    Codes cleaned and annotated: August, 2026 (XF)

    NOTE: please ignore functions labelled with [NOT USED IN THE LTM STUDY]
*/

#ifndef EVOLVETUMOUR_HPP_INCLUDED
#define EVOLVETUMOUR_HPP_INCLUDED

#include "initTumour.hpp"

// one iteration of simulation
void oneIterOverdamp(vector<Group> &vecGroup, vector<Node> &vecNode, VoxMap &voxMap,
                     vector<double> &randNormNum0, vector<double> &randNormNumLife, vector<double> &randNum01,
                     vector<unordered_map<string, double>> &mitosisToWrite,
                     vector<Node> &vecNodeBd);

void updateCoord(vector<Node> &vecNode, vector<Group> &vecGroup, vector<double> &randNormNum0);
void updateVeloc(vector<Group> &vecGroup, vector<Node> &vecNode, VoxMap &voxMap,
                 vector<double> &randNormNum0,
                 vector<Node> &vecNodeBd);

void mitosis(vector<Group> &vecGroup, vector<Node> &vecNode, vector<double> &randNormNumLife, vector<double> &randNum01,
             vector<unordered_map<string, double>> &mitosisToWrite);

void induceSubclones(vector<Group> &vecGroup, vector<Node> &vecNode, vector<double> &randNum01);

// [NOT USED IN THE LTM STUDY] following functions
void emergeSubclones(vector<Group> &vecGroup, vector<Node> &vecNode, vector<int> &vecNodeIdsDivide);
void recordMutations(unordered_map<int, unordered_map<string, int>> umap_mutations);

// VoxMap algorithm
void updateVoxMap(VoxMap &voxMap, vector<Node> &vecNode);

double dVecDist(dVec &p1, dVec &p2);

#endif
