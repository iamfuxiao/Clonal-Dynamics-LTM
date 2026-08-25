/*
    File: initTumour.hpp
    Model: particleCell -- lineage tracing model (LTM)
    Created: 10 May, 2018 (XF)
    Codes cleaned and annotated: August, 2026 (XF)

    [1] In the study Bailey_Bhargava_Fu_etal, key parameters explored in the Lineage Tracing Model (LTM)

    > GAMMA_VISCOUS - which reflects the environmental drag term: values 100 (high) and 30 (low)
    > BROWNIAN_DIFF_SCALE - which reflects a cell's intrinsic random motility term: values 0 or 1 (low), 1.6e5 (intermediate), 6e5 (high)
    > MITOSIS_T_DOUBLE_AVG - which reflects a cell's intrinsic doubling time: default value 15

    Expectations:
        lower GAMMA_VISCOUS and higher BROWNIAN_DIFF_SCALE
            attenuates growth arrest and causes more subclone fragmentation and inter-subclone mixing
        larger MITOSIS_T_DOUBLE_AVG
            also attenuates growth arrest

    [2] Other parameters are not explored/discussed extensively in the study Bailey_Bhargava_Fu_etal
        Readers and users may find these also interesting to explore using the model

        for example,
    > MITOSIS_CIP_THRP - which reflects a cell's sensitivity to compression-induced proliferative restriction
        comment: the smaller the absolute value, the thinner the active proliferative layer in the simulated colony and the rougher the colony frontier
    > NUM_HOUR_INDUCE_SUBCL & PERC_CELL_INDUCE_SUBCL - which control the timing and extent of subclone labelling
        comment: earlier timing of labelling & lower fraction of cells being labelled resulted in fewer subclones


    NOTE: please ignore parameters and functions labelled with [NOT USED IN THE LTM STUDY]
*/

#ifndef INITTUMOUR_HPP_INCLUDED
#define INITTUMOUR_HPP_INCLUDED

#include <stdio.h>
// #include <math.h>
#include <iostream>
#include <fstream>

#include <cstring>
#include <string>
#include <sstream>
#include <iomanip>

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <cmath>

#include <random>
#include <chrono>

#include <omp.h>
#include <unistd.h> // for getpid

using namespace std;

typedef struct
{
    double x, y, z;
} dVec;
typedef struct
{
    double x, y;
} dPair;
typedef struct
{
    int i, j, k;
} iVec;
typedef struct
{
    int i, j;
} iPair;

/* ~~~~~~~~~~~ Constants ~~~~~~~~~~ */
const double PI = 3.1415926535897;
const double K_BOLTZMANN = 1.38E-2; // (pN*nm/K)
const double TEMPERATURE = 310;     // (K)  body temperature
const string SIM_DIM = "2D";        // "2D" or "3D"
const pid_t PROC_ID = getpid();
const int DIM_VOXMAP = 16; // n x n x n
const int NUM_THREADS = 2;

/* ~~~~~~~~~~~ Optional Module Switchs ~~~~~~~~~~ */
const bool BROWNIAN_ON = true;    // brownian motion / stochastic force
const bool BOUNDARY_ON = false;   // [NOT USED IN THE LTM STUDY] boundary elements
const bool THERAPY_ON = false;    // [NOT USED IN THE LTM STUDY] virtual chemotherapy
const bool PROPULSION_ON = false; // [NOT USED IN THE LTM STUDY] propulsive forces
const bool APOPTOSIS_ON = false;  // [NOT USED IN THE LTM STUDY] apoptosis

/* ~~~~~~~~~~~ Parameters ~~~~~~~~~~ */
const double DT = 0.01;   // hour
const double T = 24 * 14; // hour; 24*10 is used for the demo; value of 24*14 or 24*21 used in the study
const int FREQ_PRINT = (int)(4 / DT);
const int FREQ_WRITE = (int)(4 / DT);
const int FREQ_WRITE_DYNAMICS = (int)(4 / DT);
const int FREQ_WRITE_MUTATIONS = (int)(24 * 21 / DT); // [NOT USED IN THE LTM STUDY] original value 0.1
const int FREQ_WRITE_MITOSIS = (int)(24 / DT);

const double CELL_RADIUS = 15;                         // um
const double TUMOUR_RADIUS_MAX = 4096;                 // um; stopping condition using a size threshold
const vector<int> VEC_TUMOUR_RADIUS_WRITE_MUTATIONS{}; // [NOT USED IN THE LTM STUDY]

// Initial configuration
const string INIT_SPECIAL_PATTERN = "none"; // "none", "three_line"
const int INIT_NUM_GROUP = 1;
const int INIT_NUM_NODE_PER_GROUP = 1; // 1 if special pattern is none; 60 for three_line
const double INIT_TUMOUR_RADIUS = cbrt(INIT_NUM_GROUP * INIT_NUM_NODE_PER_GROUP) * CELL_RADIUS;
const double INIT_TUMOUR_RADIUS_2D = sqrt(INIT_NUM_GROUP * INIT_NUM_NODE_PER_GROUP) * CELL_RADIUS;

// Brownian motion

// -- [ Key parameter studied in Bailey_Bhargava_Fu_etal ]
const double GAMMA_VISCOUS = 1E2; // drag coefficient
// -- [ Key parameter studied in Bailey_Bhargava_Fu_etal ]
const double BROWNIAN_DIFF_SCALE = 0;
const double BROWNIAN_DIFF_COEF_SQRT = sqrt(BROWNIAN_DIFF_SCALE * 2 * K_BOLTZMANN * TEMPERATURE * 3600 * 3600 / GAMMA_VISCOUS * 1E-9); // sqrt(um*um/h)

// Adhesion & Repulsion
const double ADHESION_LENGTH_ACT = 2 * CELL_RADIUS;   // um
const double ADHESION_LENGTH_EQU = 1.5 * CELL_RADIUS; // um
const double ADHESION_AREA_MIN = 1;                   // um^2
// ... LJ potential [NOT USED] need to debug huge repulsion following mitosis
const double ADHESION_STRENGTH_LJ = 0; // [NOT USED IN THE LTM STUDY] pN*um
// ... Harmonic potential
const double ADHESION_STRENGTH = 1000;  // kg/h^2
const double REPULSION_STRENGTH = 1000; // kg/h^2
const double REPULSION_LENGTH_ACT = 2 * CELL_RADIUS;

// Proliferation & Mitosis
const double MITOSIS_T_DOUBLE_AVG = 15;                  // h; if it's erlang, mean=ku
const double MITOSIS_T_DOUBLE_STD = 2;                   // h; if it's erlang, var=ku^2
const string MITOSIS_T_DOUBLE_DIST = "Erlang";           // "Erlang", "Gaussian"
const string MITOSIS_CIP_TYPE = "Mechanics";             // "Mechanics", "Topology", "NoFreeEdge"
const string MITOSIS_ENP_TYPE = "none";                  // [NOT USED IN THE LTM STUDY] "Inherit", "Mechanics", "Geometry"
const int MITOSIS_CIP_THR_3D = 12;                       // [NOT USED IN THE LTM STUDY]
const int MITOSIS_CIP_THR_2D = 6;                        // [NOT USED IN THE LTM STUDY] default: 6 for 2D, 12 for 3D; number of surrounding cells needs to be smaller than or equal to this number to permit mitosis
const double MITOSIS_CIP_THRP = -17 / CELL_RADIUS * 15.; //
const int PERIPHERY_MAX_NCB = 5;                         // [NOT USED IN THE LTM STUDY]

// Subclone evolution
const string TYPE_SUBCL = "induce"; // "induce" - lineage tracing model (LTM), "emerge" - random mutation model (RMM), "none"
// .. Subclone induction (chemical exposure)
const int NUM_SUBCL = 10;                      //
const string TYPE_INDUCE_SUBCL = "stochastic"; // "stochastic"
const int NUM_CELL_INDUCE_SUBCL = 200;         // [NOT USED IN THE LTM STUDY] [x] size of clone to induce subclone
const int NUM_CELL_INDUCE_SUBCL_DET = 1000;    // [NOT USED IN THE LTM STUDY] [x] size of clone to induce subclone in "deter..." type
const int NUM_HOUR_INDUCE_SUBCL = 168;         // hour; value 120 used in demo; value 168 used in Bailey_Bhargava_Fu_etal
const double PERC_CELL_INDUCE_SUBCL = 0.1;     // the fraction of cells incuded as subclones; used for stochastic sampling
// .. Subclone emergence (mutation)
const string TYPE_EMERGE_SUBCL = "random";
const double RATE_MUT_PER_DIV = 5E-3;                             // [NOT USED IN THE LTM STUDY]
const double PROB_MUT_DRIVER = 5E-4;                              // [NOT USED IN THE LTM STUDY]
const int NUM_MUT_POOL = 10000;                                   // [NOT USED IN THE LTM STUDY]
const int NUM_MUT_PER_DIV = 5;                                    // [NOT USED IN THE LTM STUDY]
const int NUM_MUT_DRIVER = (int)(PROB_MUT_DRIVER * NUM_MUT_POOL); // [NOT USED IN THE LTM STUDY]

// [NOT USED IN THE LTM STUDY]  Brownian rotation + active propulsion
const double BROWNIAN_DIFF_COEF_SQRT_PROP_ANGLE = 0.01; // sqrt(deg*deg/h)
const double PROPULSION_SIZE = 100;                     // pN

// [NOT USED IN THE LTM STUDY]  Chemotherapy
const string TYPE_DRUG = "cytotoxic";        // "cytotoxic"
const string TYPE_DRUG_SCHEDULE = "uniform"; // "uniform", "per_num_cell"
const double DRUG_DOSAGE = 1;                // 1 is max dosage
const double DRUG_PERIOD = 48;               // hours
const double DRUG_PER_NUM_CELL = 10000;      //
const int NUM_HOUR_INTRO_DRUG = DRUG_PERIOD; // hours post clonal induction
// const int NUM_HOUR_EFFECT_DRUG      =   0;  // [x] not used yet

// [NOT USED IN THE LTM STUDY]  Apoptosis
const double APOPTOSIS_RATE = 1e-3; // 1/h;

// [OPTIONAL] Group-specific Scaling coefficients
const string TYPE_DRIVER_MUTATION = "MITOSIS"; // "MITOSIS", "MITOSIS_PTHR", "ADHESION", "FRICTION", "SENSITIVITY_DRUG"
const double GROUP_SCALE_VIS = 1;              // scale GAMMA_VISCOUS
const double GROUP_SCALE_ADH = 1;              // scale ADHESION_STRENGTH
const double GROUP_SCALE_MIT = 1;              // scale MITOSIS_T_DOUBLE
const double GROUP_SCALE_SEN = 1;              // scale drug sensitivity 1 or 0

// [NOT USED IN THE LTM STUDY]  Boundary effects
const string TYPE_BOUNDARY = "radcir";                       // "random", "vertical", "radial", "circum"
const int INIT_NUM_GROUPBD = 200;                            // number of groupGd
const int INIT_NUM_NODE_PER_GROUPBD = 25;                    // number of nodes per groupBd
const double POS_BOUNDARY_LINE_X = 200;                      // location of a line-shaped boundary; only used if is single "vertical" boundary
const double REPULSION_STRENGTH_BD = 2 * REPULSION_STRENGTH; // assumpting twice the cell-cell repulsion

/* ~~~~~~~~~~~ Classes ~~~~~~~~~~ */
/* -------------- Node (representing Cell) --------------- */
class Node
{
    bool isAlive;         // isAlive is a flag to indicate whether the cell is dead (apoptotic or necrotic)
    int id, groupId, gen; // gen refers to generation; founder cell has generation 0
    int lineageId;        // every clonal founder cell has a unique lineageId
    dVec coord, coord_prev, veloc;
    double propu_angle;
    double time_double, time_double_intrinsic;
    vector<int> vecNbNodeIds;
    vector<double> vecNbNodeDists;

    vector<int> vecMutIds; // record a list of mutations accumulated in the cell

    double pressure, radius;
    vector<dPair> stress2d;
    vector<double> vecNbCommAreas;
    vector<double> vecNbPressures;
    unordered_map<string, int> umapCellDensity; // string indicates range; int indicates count

public:
    Node();
    Node(bool isAlive, int id, int lineageId, int groupId);
    Node(bool isAlive, int id, int lineageId, int groupId, int gen,
         dVec coord, dVec coord_prev, dVec veloc, double propu_angle,
         double time_double, double time_double_intrinsic,
         vector<int> vecNbNodeIds, vector<double> vecNbNodeDists,
         vector<int> vecMutIds,
         double pressure, double radius, vector<dPair> stress2d,
         vector<double> vecNbCommAreas, vector<double> vecNbPressures,
         unordered_map<string, int> umapCellDensity);
    ~Node();

    // setters
    void set_isAlive(bool isAlive) { this->isAlive = isAlive; }
    void set_id(int id) { this->id = id; }
    void set_lineageId(int lineageId) { this->lineageId = lineageId; }
    void set_groupId(int groupId) { this->groupId = groupId; }
    void set_gen(int gen) { this->gen = gen; }
    void set_coord(dVec coord) { this->coord = coord; }
    void set_coord_prev(dVec coord_prev) { this->coord_prev = coord_prev; }
    void set_veloc(dVec veloc) { this->veloc = veloc; }
    void set_propu_angle(double propu_angle) { this->propu_angle = propu_angle; }
    void set_time_double(double time_double) { this->time_double = time_double; }
    void set_time_double_intrinsic(double time_double_intrinsic) { this->time_double_intrinsic = time_double_intrinsic; }
    void set_vecNbNodeIds(vector<int> vecNbNodeIds) { this->vecNbNodeIds = vecNbNodeIds; }
    void set_vecNbNodeDists(vector<double> vecNbNodeDists) { this->vecNbNodeDists = vecNbNodeDists; }
    void set_vecMutIds(vector<int> vecMutIds) { this->vecMutIds = vecMutIds; }
    void set_pressure(double pressure) { this->pressure = pressure; }
    void set_radius(double radius) { this->radius = radius; }
    void set_stress2d(vector<dPair> stress2d) { this->stress2d = stress2d; }
    void set_vecNbCommAreas(vector<double> vecNbCommAreas) { this->vecNbCommAreas = vecNbCommAreas; }
    void set_vecNbPressures(vector<double> vecNbPressures) { this->vecNbPressures = vecNbPressures; }
    void set_umapCellDensity(unordered_map<string, int> umapCellDensity) { this->umapCellDensity = umapCellDensity; }

    // getters
    bool get_isAlive() const { return this->isAlive; }
    int get_id() const { return this->id; }
    int get_lineageId() const { return this->lineageId; }
    int get_groupId() const { return this->groupId; }
    int get_gen() const { return this->gen; }
    dVec get_coord() const { return this->coord; }
    dVec get_coord_prev() const { return this->coord_prev; }
    dVec get_veloc() const { return this->veloc; }
    double get_propu_angle() const { return this->propu_angle; }
    double get_time_double() const { return this->time_double; }
    double get_time_double_intrinsic() const { return this->time_double_intrinsic; }
    vector<int> get_vecNbNodeIds() const { return this->vecNbNodeIds; }
    vector<double> get_vecNbNodeDists() const { return this->vecNbNodeDists; }
    vector<int> get_vecMutIds() const { return this->vecMutIds; }
    double get_pressure() const { return this->pressure; }
    double get_radius() const { return this->radius; }
    vector<dPair> get_stress2d() const { return this->stress2d; }
    vector<double> get_vecNbCommAreas() const { return this->vecNbCommAreas; }
    vector<double> get_vecNbPressures() const { return this->vecNbPressures; }
    unordered_map<string, int> get_umapCellDensity() const { return this->umapCellDensity; }
};

/* -------------- Group (representing Subclone) ------------- */
class Group
{
    int id;
    vector<int> vecNodeIds;
    double scale_vis, scale_adh, scale_mit, scale_sen; // this is scaling coefficient on the global parameter value

public:
    Group();
    Group(int id, vector<int> vecNodeIds, double scale_vis, double scale_adh, double scale_mit, double scale_sen);
    ~Group();

    // setters
    void set_id(int id) { this->id = id; }
    void set_vecNodeIds(vector<int> vecNodeIds) { this->vecNodeIds = vecNodeIds; }
    void set_scale_vis(double scale_vis) { this->scale_vis = scale_vis; }
    void set_scale_adh(double scale_adh) { this->scale_adh = scale_adh; }
    void set_scale_mit(double scale_mit) { this->scale_mit = scale_mit; }
    void set_scale_sen(double scale_sen) { this->scale_sen = scale_sen; }

    // getters
    int get_id() const { return this->id; }
    vector<int> get_vecNodeIds() const { return this->vecNodeIds; }
    double get_scale_vis() const { return this->scale_vis; }
    double get_scale_adh() const { return this->scale_adh; }
    double get_scale_mit() const { return this->scale_mit; }
    double get_scale_sen() const { return this->scale_sen; }
};

/* -------------- VoxMap ----------------- */
class VoxMap
{
    iVec dim;
    unordered_map<int, vector<int>> mapBeadIds;     // contains Bead ids in VoxBox; key is collapsed from 3d indexing
    unordered_map<int, vector<int>> mapBeadIdsHalo; // contains Bead ids in halo of VoxBox;

public:
    VoxMap();
    VoxMap(iVec dim, unordered_map<int, vector<int>> mapBeadIds, unordered_map<int, vector<int>> mapBeadIdsHalo);
    ~VoxMap();

    // setters
    void set_dim(iVec dim) { this->dim = dim; }
    void set_mapBeadIds(unordered_map<int, vector<int>> mapBeadIds) { this->mapBeadIds = mapBeadIds; }
    void set_mapBeadIdsHalo(unordered_map<int, vector<int>> mapBeadIdsHalo) { this->mapBeadIdsHalo = mapBeadIdsHalo; }

    // getters
    iVec get_dim() const { return this->dim; }
    unordered_map<int, vector<int>> get_mapBeadIds() const { return this->mapBeadIds; }
    unordered_map<int, vector<int>> get_mapBeadIdsHalo() const { return this->mapBeadIdsHalo; }
};

/* ~~~~~~~~~~~ Functions ~~~~~~~~~~ */
void initGroup(vector<Group> &vecGroup, vector<Node> &vecNode);
void initGroupBd(vector<Group> &vecGroupBd, vector<Node> &vecNodeBd);
void initDynamics(vector<Group> &vecGroup, vector<Node> &vecNode, vector<double> &randNum01, vector<double> &randNormNumLife);
void initDynamicsSpecialPattern(vector<Group> &vecGroup, vector<Node> &vecNode, vector<double> &randNormNumLife, string typePattern);
void createLineOfCoords(vector<dVec> &vecLineCoord, int ncell, dVec orientation, double stepsize);
void initDynamicsBd(vector<Group> &vecGroupBd, vector<Node> &vecNodeBd, vector<double> &randNum01);

void currCloneCenterPeriphery(vector<Node> &vecNode, vector<int> &vecNodePeriphery, dVec &cloneCOM, bool flagCheckPeriphery);
void sortClonePeriphery(vector<Node> &vecNode, vector<int> &vecNodePeriphery);
void fillGapPeriphery(vector<Node> &vecNode, vector<int> &vecNodePeriphery);
void calcCurvature(vector<Node> &vecNode, vector<int> &vecNodePeriphery, unordered_map<int, double> &mapNodeCurvature);

double printGroupInfo(vector<Group> &vecGroup, vector<Node> &vecNode, double t_now, int seed);
void writeGroupInfo(vector<Group> &vecGroup, vector<Node> &vecNode, double t_now, vector<Node> &vecNodeBd);
void writeNodeDynamics(vector<Node> &vecNode, double t_now, int seed, bool flagWriteMut);
void writeNodeMitosis(vector<unordered_map<string, double>> &mitosisToWrite, double t_now);

double dVecDist(dVec &p1, dVec &p2);

#endif
