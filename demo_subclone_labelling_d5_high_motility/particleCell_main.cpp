/*
    File: particleCell_main.cpp
    Model: particleCell -- lineage tracing model (LTM)
    Created: 10 May, 2018 (XF)
    Codes cleaned and annotated: August, 2026 (XF)

    NOTE: please ignore functions labelled with [NOT USED IN THE LTM STUDY]

*/

#include "initTumour.hpp"
#include "evolveTumour.hpp"

int main(int argc, char **argv)
{

    // set up random number generator
    // int seed = chrono::system_clock::now().time_since_epoch().count();
    int seed = 123;
    mt19937 generator(seed);
    uniform_real_distribution<double> uniform01(0.0, 1.0);
    normal_distribution<double> normal0(0.0, sqrt(DT));
    normal_distribution<double> normalT(MITOSIS_T_DOUBLE_AVG, MITOSIS_T_DOUBLE_STD);

    // initialize voxMap
    const int dimX = DIM_VOXMAP;
    iVec dim = {dimX, dimX, dimX};
    if (SIM_DIM == "2D")
        iVec dim = {dimX, dimX, 1};
    unordered_map<int, vector<int>> mapBeadIds;
    unordered_map<int, vector<int>> mapBeadIdsHalo;
    VoxMap voxMap = VoxMap(dim, mapBeadIds, mapBeadIdsHalo);

    // initialize vectors for writing to file
    vector<unordered_map<string, double>> mitosisToWrite;

    // initialize vectors
    vector<Group> vecGroup;   // collect all Groups (i.e., Subclones)
    vector<Node> vecNode;     // collect all Nodes (i.e., Cells)
    vector<Group> vecGroupBd; // collect all Groups (i.e., different Boundaries: stroma, bronchiole, vasculature, ... Or simply container wall)
    vector<Node> vecNodeBd;   // collect all Nodes (i.e., components of the GroupBd)

    // initialize groups
    initGroup(vecGroup, vecNode);

    // initialize node coordinates & velocties
    vector<double> randNum01, randNormNumLife;
    for (int i = 0; i < vecNode.size() * 3 + 1; i++)
    {
        randNum01.push_back(uniform01(generator));
        if (i < vecNode.size() + 1)
        {
            if (MITOSIS_T_DOUBLE_DIST == "Gaussian")
                randNormNumLife.push_back(normalT(generator));
            if (MITOSIS_T_DOUBLE_DIST == "Erlang") // default
            {
                double erl, factor = 1;
                const double u_erl = MITOSIS_T_DOUBLE_STD * MITOSIS_T_DOUBLE_STD / MITOSIS_T_DOUBLE_AVG;
                const int k_erl = (int)(MITOSIS_T_DOUBLE_AVG / u_erl);
                for (int kk = 0; kk < k_erl; kk++)
                    factor *= uniform01(generator);
                erl = -u_erl * log(factor);
                randNormNumLife.push_back(erl);
            }
        }
    }

    // default
    if (INIT_SPECIAL_PATTERN == "none")
        initDynamics(vecGroup, vecNode, randNum01, randNormNumLife);

    // [NOT USED IN THE LTM STUDY]
    if (INIT_SPECIAL_PATTERN == "three_line")
        initDynamicsSpecialPattern(vecGroup, vecNode, randNormNumLife, INIT_SPECIAL_PATTERN);

    // [NOT USED IN THE LTM STUDY] initialize boundary
    // ... "Boundary" or "Bd" related codes were developed to simulate cancer cell interactions with rigit structures;
    // ... these are not studied in Bailey_Bhargava_Fu_etal
    if (BOUNDARY_ON)
    {
        // initialize groupBd
        initGroupBd(vecGroupBd, vecNodeBd);
        // initialize nodeBd coordinates
        randNum01.clear();
        for (int i = 0; i < vecNodeBd.size() * 12 + 1; i++)
            randNum01.push_back(uniform01(generator));
        initDynamicsBd(vecGroupBd, vecNodeBd, randNum01);
    }

    // SIMULATION starts
    double t_now = 0;
    long nIter = 0;
    vector<double> vecTumourSize;
    while (t_now <= T)
    {
        // one iteration of simulation
        vector<double> randNormNum0;
        randNormNumLife.clear();
        randNum01.clear();
        for (int i = 0; i < vecNode.size() * 7; i++)
        {
            randNormNum0.push_back(normal0(generator));
            if (i < vecNode.size() + 1)
            {
                randNum01.push_back(uniform01(generator));
                if (MITOSIS_T_DOUBLE_DIST == "Gaussian")
                    randNormNumLife.push_back(normalT(generator));
                if (MITOSIS_T_DOUBLE_DIST == "Erlang") // default
                {
                    double erl, factor = 1;
                    const double u_erl = MITOSIS_T_DOUBLE_STD * MITOSIS_T_DOUBLE_STD / MITOSIS_T_DOUBLE_AVG;
                    const int k_erl = (int)(MITOSIS_T_DOUBLE_AVG / u_erl);
                    for (int kk = 0; kk < k_erl; kk++)
                        factor *= uniform01(generator);
                    erl = -u_erl * log(factor);
                    randNormNumLife.push_back(erl);
                }
            }
        }
        int nNode_old = vecNode.size();
        oneIterOverdamp(vecGroup, vecNode, voxMap, randNormNum0, randNormNumLife, randNum01, mitosisToWrite,
                        vecNodeBd);
        int nNode_new = vecNode.size();

        // introduce subclones in the Lineage Tracing Model (LTM)
        if (TYPE_SUBCL == "induce")
        {
            // induce subclones
            int sizeColonyInduce = NUM_CELL_INDUCE_SUBCL;
            if (TYPE_INDUCE_SUBCL == "deterministic")
                sizeColonyInduce = NUM_CELL_INDUCE_SUBCL_DET;
            // at a defined time
            if (t_now >= NUM_HOUR_INDUCE_SUBCL && t_now - DT < NUM_HOUR_INDUCE_SUBCL)
            {
                cout << "\n Inducing Subclones at t_now = " << t_now << " ! \n"
                     << endl;
                randNum01.clear();
                for (int i = 0; i < nNode_new * 10; i++)
                    randNum01.push_back(uniform01(generator));
                induceSubclones(vecGroup, vecNode, randNum01);
            }
        }

        // print information
        if (nIter % FREQ_PRINT == 0)
        {
            stringstream PROC_ID_SS;
            PROC_ID_SS << PROC_ID;
            cout << "\n ----- pid = " << PROC_ID_SS.str() << " ----- t_now = " << t_now << " -----" << endl;

            double tumourSize;
            tumourSize = printGroupInfo(vecGroup, vecNode, t_now, seed);
            vecTumourSize.push_back(tumourSize);
        }

        // [older way of saving outputs in PDB format for visualisation using PyMol] write coordinate to PDB file
        // if (nIter % FREQ_WRITE == 0)
        // {
        //     cout << "... writing group information ..." << endl;
        //     writeGroupInfo(vecGroup, vecNode, t_now, vecNodeBd);
        //     cout << "done !" << endl;
        // }

        // write coordinate & velocity to file
        if (nIter % FREQ_WRITE_DYNAMICS == 0)
        {
            bool flagWriteMut = false;
            bool flagReachTumourSize = false;
            for (vector<int>::const_iterator it_s = VEC_TUMOUR_RADIUS_WRITE_MUTATIONS.begin(); it_s != VEC_TUMOUR_RADIUS_WRITE_MUTATIONS.end(); it_s++)
            {
                int rad0 = *it_s;
                if (vecTumourSize[vecTumourSize.size() - 2] / 2 < rad0 &&
                    vecTumourSize[vecTumourSize.size() - 1] / 2 >= rad0)
                {
                    flagReachTumourSize = true;
                    break;
                }
            }

            if (nIter % FREQ_WRITE_MUTATIONS == 0 || flagReachTumourSize == true)
                flagWriteMut = true;

            // if (nIter % FREQ_WRITE_MUTATIONS == 0 ||
            //     (vecTumourSize.size() > 2 &&
            //     vecTumourSize[vecTumourSize.size()-2]/2 < TUMOUR_RADIUS_WRITE_MUTATIONS &&
            //     vecTumourSize[vecTumourSize.size()-1]/2 >= TUMOUR_RADIUS_WRITE_MUTATIONS) )
            //     flagWriteMut = true;

            // cout << "... writing node dynamics information ..." << endl;
            writeNodeDynamics(vecNode, t_now, seed, flagWriteMut);
            // cout << "done !" << endl;
        }

        // [NOT USED IN THE LTM STUDY] write mitosis information to file
        // if (nIter % FREQ_WRITE_MITOSIS == 0)
        // {
        //     cout << "... writing node mitosis information ..." << endl;
        //     writeNodeMitosis(mitosisToWrite, t_now);
        //     mitosisToWrite.clear();
        //     cout << "done !" << endl;
        // }

        t_now += DT;
        nIter++;
    }

    return 0;
}
