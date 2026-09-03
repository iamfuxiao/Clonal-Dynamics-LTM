/*
    File: initTumour.cpp
    Model: particleCell -- lineage tracing model (LTM)
    Created: 10 May, 2018 (XF)
    Codes cleaned and annotated: August, 2026 (XF)

    NOTE: please ignore functions labelled with [NOT USED IN THE LTM STUDY]
*/

#include "initTumour.hpp"

/* ~~~~~~~~~~~ Classes ~~~~~~~~~~ */
/* ------------ Node ------------- */
Node::Node() {}
Node::Node(bool isAlive, int id, int lineageId, int groupId)
{
    this->isAlive = isAlive;
    this->id = id;
    this->lineageId = lineageId;
    this->groupId = groupId;
}
Node::Node(bool isAlive, int id, int lineageId, int groupId, int gen,
           dVec coord, dVec coord_prev, dVec veloc, double propu_angle,
           double time_double, double time_double_intrinsic,
           vector<int> vecNbNodeIds, vector<double> vecNbNodeDists,
           vector<int> vecMutIds,
           double pressure, double radius, vector<dPair> stress2d,
           vector<double> vecNbCommAreas, vector<double> vecNbPressures,
           unordered_map<string, int> umapCellDensity)
{
    this->isAlive = isAlive;
    this->id = id;
    this->lineageId = lineageId;
    this->groupId = groupId;

    this->gen = gen;
    this->coord = coord;
    this->coord_prev = coord_prev;
    this->veloc = veloc;
    this->propu_angle = propu_angle;

    this->time_double = time_double;
    this->time_double_intrinsic = time_double_intrinsic;
    this->vecNbNodeIds = vecNbNodeIds;
    this->vecNbNodeDists = vecNbNodeDists;

    this->vecMutIds = vecMutIds;

    this->pressure = pressure;
    this->radius = radius;
    this->stress2d = stress2d;
    this->vecNbCommAreas = vecNbCommAreas;
    this->vecNbPressures = vecNbPressures;
    this->umapCellDensity = umapCellDensity;
}
Node::~Node() {}

/* ------------- Group -------------- */
Group::Group() {}
Group::Group(int id, vector<int> vecNodeIds, double scale_vis, double scale_adh, double scale_mit, double scale_sen)
{
    this->id = id;
    this->vecNodeIds = vecNodeIds;
    this->scale_vis = scale_vis;
    this->scale_adh = scale_adh;
    this->scale_mit = scale_mit;
    this->scale_sen = scale_sen;
}
Group::~Group() {}

/* ------------- VoxMap ---------------- */
VoxMap::VoxMap() {}
VoxMap::VoxMap(iVec dim, unordered_map<int, vector<int>> mapBeadIds, unordered_map<int, vector<int>> mapBeadIdsHalo)
{
    this->dim = dim;
    this->mapBeadIds = mapBeadIds;
    this->mapBeadIdsHalo = mapBeadIdsHalo;
}
VoxMap::~VoxMap() {}

/* ~~~~~~~~~~~ Functions ~~~~~~~~~~ */
// function initGroup() is a one-time function
void initGroup(vector<Group> &vecGroup, vector<Node> &vecNode)
{
    bool isAlive = true;
    int currGroupId = vecGroup.size() + 1;
    int currNodeId = vecNode.size();
    int currLineageId = 0;
    double scale_vis = 1., scale_adh = 1., scale_mit = 1., scale_sen = 1.;

    for (int iG = 0; iG < INIT_NUM_GROUP; iG++)
    {
        vector<int> vecNodeIds;
        for (int iN = 0; iN < INIT_NUM_NODE_PER_GROUP; iN++)
        {
            Node node0 = Node(isAlive, currNodeId, currLineageId, currGroupId);
            vecNode.push_back(node0);
            vecNodeIds.push_back(currNodeId);
            currNodeId++;
        }
        Group group0 = Group(currGroupId, vecNodeIds, scale_vis, scale_adh, scale_mit, scale_sen);
        vecGroup.push_back(group0);

        cout << "... init Group ID = " << currGroupId << " ..." << endl;
        cout << "with scale_vis = " << vecGroup[currGroupId - 1].get_scale_vis() << endl;
        cout << "     scale_adh = " << vecGroup[currGroupId - 1].get_scale_adh() << endl;
        cout << "     scale_mit = " << vecGroup[currGroupId - 1].get_scale_mit() << endl;
        cout << "     scale_sen = " << vecGroup[currGroupId - 1].get_scale_sen() << endl;

        currGroupId++;
    }
}
// [NOT USED IN THE LTM STUDY] function initGroupBd() is a one-time function that creates initial configuration of boundary
void initGroupBd(vector<Group> &vecGroupBd, vector<Node> &vecNodeBd)
{
    bool isAlive = true;
    int currGroupId = vecGroupBd.size() + 1;
    int currNodeId = vecNodeBd.size();
    int currLineageId = -1;
    double scale_vis = 1., scale_adh = 1., scale_mit = 1., scale_sen = 1.;

    // const int INIT_NUM_NODE_PER_GROUPBD = 100;

    // for (int iG = 0; iG < INIT_NUM_GROUP; iG ++)
    for (int iG = 0; iG < INIT_NUM_GROUPBD; iG++)
    {
        vector<int> vecNodeIds;
        for (int iN = 0; iN < INIT_NUM_NODE_PER_GROUPBD; iN++)
        {
            Node node0 = Node(isAlive, currNodeId, currLineageId, currGroupId);
            vecNodeBd.push_back(node0);
            vecNodeIds.push_back(currNodeId);
            currNodeId++;
        }
        Group group0 = Group(currGroupId, vecNodeIds, scale_vis, scale_adh, scale_mit, scale_sen);
        vecGroupBd.push_back(group0);

        cout << "... init GroupBd ID = " << currGroupId << " ..." << endl;
        currGroupId++;
    }
}

// function initDynamics() is a one-time function
void initDynamics(vector<Group> &vecGroup, vector<Node> &vecNode, vector<double> &randNum01, vector<double> &randNormNumLife)
{
    double rInit = INIT_TUMOUR_RADIUS;
    double rInit2D = INIT_TUMOUR_RADIUS_2D;
    double r, theta, phi;
    int k = 0, kk = 0;

    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it++)
    {
        // coordinate
        dVec coord;
        if (SIM_DIM == "3D")
        {
            r = rInit * cbrt(randNum01[k++]);
            theta = acos(2. * randNum01[k++] - 1);
            phi = 2. * PI * randNum01[k++];
            coord = {r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta)};
        }
        if (SIM_DIM == "2D")
        {
            r = rInit2D * sqrt(randNum01[k++]);
            phi = 2. * PI * randNum01[k++];
            coord = {r * cos(phi), r * sin(phi), 0.};
        }
        it->set_coord(coord);
        it->set_coord_prev(coord);

        // velocity
        it->set_veloc({0., 0., 0.});

        // [NOT USED IN THE LTM STUDY] propulsion angle
        // double propu_angle = rand() / (float)RAND_MAX * 2 * PI;
        double propu_angle = 0;
        it->set_propu_angle(propu_angle);
        // cout << propu_angle << endl;

        // time_double
        double time_double = randNormNumLife[kk++];
        it->set_time_double(time_double);
        it->set_time_double_intrinsic(time_double);

        // generation
        it->set_gen(0);

        // accumulate mutations
        vector<int> vecMutIds;
        it->set_vecMutIds(vecMutIds);

        // mechanics
        it->set_pressure(0.);
        it->set_radius(CELL_RADIUS);
    }
}
// [NOT USED IN THE LTM STUDY] function initDynamicsSpecialPattern() is a one-time function
void initDynamicsSpecialPattern(vector<Group> &vecGroup, vector<Node> &vecNode, vector<double> &randNormNumLife, string typePattern)
{
    if (typePattern == "three_line")
    {
        // (1) create three cells to form a triangle
        // (2) elongate these cells along a line
        int nCellPerLine = vecNode.size() / 3;
        vector<dVec> vecLineCoord1, vecLineCoord2, vecLineCoord3;
        // ... verticle line
        vecLineCoord1.push_back({0, sqrt(3) / 3. * 2 * CELL_RADIUS, 0});
        createLineOfCoords(vecLineCoord1, nCellPerLine, {0, 1, 0}, 2 * CELL_RADIUS);
        // ... lower left line
        vecLineCoord2.push_back({-CELL_RADIUS, -sqrt(3) / 6. * 2 * CELL_RADIUS, 0});
        createLineOfCoords(vecLineCoord2, nCellPerLine, {-sqrt(3) / 2, -0.5, 0}, 2 * CELL_RADIUS);
        // ... lower right line
        vecLineCoord3.push_back({CELL_RADIUS, -sqrt(3) / 6. * 2 * CELL_RADIUS, 0});
        createLineOfCoords(vecLineCoord3, vecNode.size() - 2 * nCellPerLine, {sqrt(3) / 2, -0.5, 0}, 2 * CELL_RADIUS);

        int kk = 0;
        for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it++)
        {
            dVec coord;
            if (it - vecNode.begin() < nCellPerLine)
                coord = vecLineCoord1[it - vecNode.begin()];
            else if (it - vecNode.begin() < 2 * nCellPerLine)
                coord = vecLineCoord2[(it - vecNode.begin()) % nCellPerLine];
            else
                coord = vecLineCoord3[(it - vecNode.begin()) % (2 * nCellPerLine)];

            it->set_coord(coord);
            it->set_coord_prev(coord);

            // velocity
            it->set_veloc({0., 0., 0.});

            // time_double
            double time_double = randNormNumLife[kk++];
            it->set_time_double(time_double);
            it->set_time_double_intrinsic(time_double);

            // generation
            it->set_gen(0);

            // mechanics
            it->set_pressure(0.);
            it->set_radius(CELL_RADIUS);
        }
    }
}
// [NOT USED IN THE LTM STUDY]
void createLineOfCoords(vector<dVec> &vecLineCoord, int ncell, dVec orientation, double stepsize)
{
    // Note: ensure that vecLineCoord has 1 coord already!
    int cnt = 1;
    while (cnt < ncell)
    {
        dVec coord0 = vecLineCoord[cnt - 1];
        dVec coord1 = {coord0.x + orientation.x * stepsize,
                       coord0.y + orientation.y * stepsize,
                       coord0.z + orientation.z * stepsize};
        vecLineCoord.push_back(coord1);
        cnt++;
    }
}
// [NOT USED IN THE LTM STUDY] function initDynamicsBd() is a one-time function that initialize coordinates of nodeBd
void initDynamicsBd(vector<Group> &vecGroupBd, vector<Node> &vecNodeBd, vector<double> &randNum01)
{
    // Assumption 1 : a straight line of nodes
    if (TYPE_BOUNDARY == "vertical")
    {
        // const double POS_BOUNDARY_LINE_X = 200;
        int nCellPerLine = vecNodeBd.size();
        double y_min = -CELL_RADIUS * nCellPerLine;
        vector<dVec> vecLineCoord;
        vecLineCoord.push_back({POS_BOUNDARY_LINE_X, y_min, 0});
        createLineOfCoords(vecLineCoord, nCellPerLine, {0, 1, 0}, 2 * CELL_RADIUS);

        int kk = 0;
        for (vector<Node>::iterator it = vecNodeBd.begin(); it != vecNodeBd.end(); it++)
        {
            dVec coord;
            if (it - vecNodeBd.begin() < nCellPerLine)
                coord = vecLineCoord[it - vecNodeBd.begin()];

            it->set_coord(coord);
            it->set_coord_prev(coord);

            // velocity
            it->set_veloc({0., 0., 0.});

            // time_double  --   these NodeBds don't divide
            double time_double = T + DT;
            it->set_time_double(time_double);
            it->set_time_double_intrinsic(time_double);

            // generation
            it->set_gen(0);

            // mechanics
            it->set_pressure(0.);
            it->set_radius(CELL_RADIUS); // assuming the same as cancer cells
        }
    }
    if (TYPE_BOUNDARY == "random" || TYPE_BOUNDARY == "radial" || TYPE_BOUNDARY == "circum" || TYPE_BOUNDARY == "radcir")
    {
        // need to set coordinates for all groupBd
        int k = 0;
        for (vector<Group>::iterator it_g = vecGroupBd.begin(); it_g != vecGroupBd.end(); it_g++)
        {
            vector<int> vecNodeIds = it_g->get_vecNodeIds();
            int nCellThisLine = vecNodeIds.size();
            double x_start, y_start, dx_ori, dy_ori;

            // randomize the starting point
            // x_start = (2.* ((double) rand() / (RAND_MAX)) - 1) * TUMOUR_RADIUS_MAX;
            // y_start = (2.* ((double) rand() / (RAND_MAX)) - 1) * TUMOUR_RADIUS_MAX;
            x_start = (2. * randNum01[k++] - 1) * TUMOUR_RADIUS_MAX;
            y_start = (2. * randNum01[k++] - 1) * TUMOUR_RADIUS_MAX;
            while (x_start * x_start + y_start * y_start > TUMOUR_RADIUS_MAX * TUMOUR_RADIUS_MAX)
            {
                // x_start = (2.* ((double) rand() / (RAND_MAX)) - 1) * TUMOUR_RADIUS_MAX;
                // y_start = (2.* ((double) rand() / (RAND_MAX)) - 1) * TUMOUR_RADIUS_MAX;
                x_start = (2. * randNum01[k++] - 1) * TUMOUR_RADIUS_MAX;
                y_start = (2. * randNum01[k++] - 1) * TUMOUR_RADIUS_MAX;
            }
            // dx_ori = 2.* ((double) rand() / (RAND_MAX)) - 1;
            // dy_ori = 2.* ((double) rand() / (RAND_MAX)) - 1;

            // This is random orientation
            dx_ori = 2. * randNum01[k++] - 1;
            dy_ori = 2. * randNum01[k++] - 1;
            while (dx_ori * dx_ori + dy_ori * dy_ori > 1)
            {
                dx_ori = 2. * randNum01[k++] - 1;
                dy_ori = 2. * randNum01[k++] - 1;
            }
            double ori_size = sqrt(dx_ori * dx_ori + dy_ori * dy_ori);
            dx_ori /= ori_size;
            dy_ori /= ori_size;
            //
            // This is radial orientation
            if (TYPE_BOUNDARY == "radial")
            {
                double theta = atan2(y_start, x_start);
                dx_ori = cos(theta);
                dy_ori = sin(theta);
            }
            //
            // This is circum orientation
            if (TYPE_BOUNDARY == "circum")
            {
                double theta = atan2(y_start, x_start);
                dy_ori = cos(theta);
                dx_ori = -sin(theta);
            }
            //
            // This is radial (left) - circum (right)
            if (TYPE_BOUNDARY == "radcir")
            {
                double theta = atan2(y_start, x_start);
                dx_ori = cos(theta); // left  : radial
                dy_ori = sin(theta);
                if (x_start > 0) // right : circum
                {
                    dy_ori = cos(theta);
                    dx_ori = -sin(theta);
                }
            }

            vector<dVec> vecLineCoord;
            vecLineCoord.push_back({x_start, y_start, 0});
            createLineOfCoords(vecLineCoord, nCellThisLine, {dx_ori, dy_ori, 0}, 1.5 * CELL_RADIUS);
            for (vector<int>::iterator it = vecNodeIds.begin(); it != vecNodeIds.end(); it++)
            {
                dVec coord;
                if (it - vecNodeIds.begin() < nCellThisLine)
                    coord = vecLineCoord[it - vecNodeIds.begin()];

                int nodeId = *it;
                vecNodeBd[nodeId].set_coord(coord);
                vecNodeBd[nodeId].set_coord_prev(coord);

                // velocity
                vecNodeBd[nodeId].set_veloc({0., 0., 0.});

                // time_double  --   these NodeBds don't divide
                double time_double = T + DT;
                vecNodeBd[nodeId].set_time_double(time_double);
                vecNodeBd[nodeId].set_time_double_intrinsic(time_double);

                // generation
                vecNodeBd[nodeId].set_gen(0);

                // mechanics
                vecNodeBd[nodeId].set_pressure(0.);
                vecNodeBd[nodeId].set_radius(CELL_RADIUS); // assuming the same as cancer cells
            }
        }
    }
}

// [NOT USED IN THE LTM STUDY] the following functions
void currCloneCenterPeriphery(vector<Node> &vecNode, vector<int> &vecNodePeriphery, dVec &cloneCOM, bool flagCheckPeriphery)
{
    double xall = 0, yall = 0, zall = 0;
    int nNode = vecNode.size();
    vector<int> vecNodePeriphery_tmp;
    for (vector<Node>::iterator it_n = vecNode.begin(); it_n != vecNode.end(); it_n++)
    {
        dVec coord = it_n->get_coord();
        xall += coord.x;
        yall += coord.y;
        zall += coord.z;

        // method 2: decide periphery by umapCellDensity
        unordered_map<string, int> umapCellDensity = it_n->get_umapCellDensity();
        if (umapCellDensity["d1"] <= PERIPHERY_MAX_NCB)
            vecNodePeriphery_tmp.push_back(it_n->get_id());
    }
    cloneCOM = {xall / nNode, yall / nNode, zall / nNode};

    if (flagCheckPeriphery) // this removes nodes that are wrongly included
    {
        // cout << "number of candidate nodes at colony periphery : " << vecNodePeriphery_tmp.size() << endl;

        // testing by check inter-cell distances within periphery
        // ... problem: including many internal cells
        if (false)
        {
            for (vector<int>::iterator it_i = vecNodePeriphery_tmp.begin(); it_i != vecNodePeriphery_tmp.end(); it_i++)
            {
                bool flagIsNodePeriphery = false;
                double dist_test = 3 * CELL_RADIUS;

                dVec coord = vecNode[*it_i].get_coord();
                // dVec coord = vecNode[*it_i].get_coord_prev();
                for (vector<int>::iterator it_j = vecNodePeriphery_tmp.begin(); it_j != vecNodePeriphery_tmp.end(); it_j++)
                {
                    if (*it_i != *it_j)
                    {
                        dVec coord2 = vecNode[*it_j].get_coord();
                        // dVec coord2 = vecNode[*it_j].get_coord_prev();
                        double dist = dVecDist(coord, coord2);
                        if (dist <= dist_test)
                        {
                            // cout << endl;
                            // cout << *it_i << " with dist = " << dist/CELL_RADIUS << endl;
                            // cout << *it_i << ": " << coord.x/CELL_RADIUS  << ", " << coord.y/CELL_RADIUS  << endl;
                            // cout << *it_j << ": " << coord2.x/CELL_RADIUS << ", " << coord2.y/CELL_RADIUS << endl;
                            flagIsNodePeriphery = true;
                            break;
                        }
                    }
                }

                if (flagIsNodePeriphery)
                    vecNodePeriphery.push_back(*it_i);
            }
        }

        // [BEST SO FAR: testing by moving cell radially
        // ... problem: miss some cells near protrusive structures
        if (true)
        {
            for (vector<int>::iterator it_i = vecNodePeriphery_tmp.begin(); it_i != vecNodePeriphery_tmp.end(); it_i++)
            {
                dVec coord = vecNode[*it_i].get_coord();
                // dVec coord = vecNode[*it_i].get_coord_prev();
                double dist = dVecDist(coord, cloneCOM), factorTest = (3 * CELL_RADIUS + dist) / dist;
                dVec vecTest = {(coord.x - cloneCOM.x) * factorTest,
                                (coord.y - cloneCOM.y) * factorTest,
                                (coord.z - cloneCOM.z) * factorTest};
                dVec pointTest = {cloneCOM.x + vecTest.x, cloneCOM.y + vecTest.y, cloneCOM.z + vecTest.z};

                bool flagIsNodePeriphery = true;
                for (vector<Node>::iterator it_n = vecNode.begin(); it_n != vecNode.end(); it_n++)
                {
                    int nId2 = it_n->get_id();
                    // if (find(vecNodePeriphery_tmp.begin(),vecNodePeriphery_tmp.end(),nId2) == vecNodePeriphery_tmp.end())
                    if (nId2 != *it_i)
                    {
                        dVec coord2 = it_n->get_coord();
                        // dVec coord2 = it_n->get_coord_prev();
                        double dist2 = dVecDist(coord2, pointTest);
                        if (dist2 < 2 * CELL_RADIUS)
                        {
                            flagIsNodePeriphery = false;
                            // cout << "... candidate ruled out : node id = " << *it_i << endl;
                            break;
                        }
                    }
                }

                if (flagIsNodePeriphery)
                    vecNodePeriphery.push_back(*it_i);
            }
        }

        // testing by moving cell towards empty edge (decided according to instantaneous force)
        // ... problem: wrong when cell has tensile forces
        if (false)
        {
            for (vector<int>::iterator it_i = vecNodePeriphery_tmp.begin(); it_i != vecNodePeriphery_tmp.end(); it_i++)
            {
                dVec coord = vecNode[*it_i].get_coord();
                dVec veloc = vecNode[*it_i].get_veloc();
                double speed = sqrt(veloc.x * veloc.x + veloc.y * veloc.y + veloc.z * veloc.z);
                double factorTest = 3 * CELL_RADIUS;
                dVec vecTest = {
                    veloc.x / speed * factorTest, veloc.y / speed * factorTest, veloc.z / speed * factorTest};
                dVec pointTest = {
                    coord.x + vecTest.x, coord.y + vecTest.y, coord.z + vecTest.z};

                bool flagIsNodePeriphery = true;
                for (vector<Node>::iterator it_n = vecNode.begin(); it_n != vecNode.end(); it_n++)
                {
                    int nId2 = it_n->get_id();
                    // if (find(vecNodePeriphery_tmp.begin(),vecNodePeriphery_tmp.end(),nId2) == vecNodePeriphery_tmp.end())
                    if (nId2 != *it_i)
                    {
                        dVec coord2 = it_n->get_coord();
                        // dVec coord2 = it_n->get_coord_prev();
                        double dist2 = dVecDist(coord2, pointTest);
                        if (dist2 < 2 * CELL_RADIUS)
                        {
                            flagIsNodePeriphery = false;
                            // cout << "... candidate ruled out : node id = " << *it_i << endl;
                            break;
                        }
                    }
                }

                if (flagIsNodePeriphery)
                    vecNodePeriphery.push_back(*it_i);
            }
        }

        // testing by moving cell towards empty edge (decided according to neighbors)
        // ...
        if (true)
        {
            for (vector<int>::iterator it_i = vecNodePeriphery_tmp.begin(); it_i != vecNodePeriphery_tmp.end(); it_i++)
            {
                dVec coord = vecNode[*it_i].get_coord();
                vector<int> vecNbNodeIds = vecNode[*it_i].get_vecNbNodeIds();
                dVec sumVec = {0, 0, 0};
                for (vector<int>::iterator it_j = vecNbNodeIds.begin(); it_j != vecNbNodeIds.end(); it_j++)
                {
                    dVec coord2 = vecNode[*it_j].get_coord();
                    dVec sumVec_tmp = sumVec;
                    sumVec = {
                        sumVec_tmp.x + coord.x - coord2.x,
                        sumVec_tmp.y + coord.y - coord2.y,
                        sumVec_tmp.z + coord.z - coord2.z};
                }
                double sumVecSize = sqrt(sumVec.x * sumVec.x + sumVec.y * sumVec.y + sumVec.z * sumVec.z);

                double factorTest = 3.5 * CELL_RADIUS;
                dVec vecTest = {
                    sumVec.x / sumVecSize * factorTest, sumVec.y / sumVecSize * factorTest, sumVec.z / sumVecSize * factorTest};
                dVec pointTest = {
                    coord.x + vecTest.x, coord.y + vecTest.y, coord.z + vecTest.z};

                bool flagIsNodePeriphery = true;
                for (vector<Node>::iterator it_n = vecNode.begin(); it_n != vecNode.end(); it_n++)
                {
                    int nId2 = it_n->get_id();
                    // if (find(vecNodePeriphery_tmp.begin(),vecNodePeriphery_tmp.end(),nId2) == vecNodePeriphery_tmp.end())
                    if (nId2 != *it_i)
                    {
                        dVec coord2 = it_n->get_coord();
                        // dVec coord2 = it_n->get_coord_prev();
                        double dist2 = dVecDist(coord2, pointTest);
                        if (dist2 < 2 * CELL_RADIUS)
                        {
                            flagIsNodePeriphery = false;
                            // cout << "... candidate ruled out : node id = " << *it_i << endl;
                            break;
                        }
                    }
                }

                if (flagIsNodePeriphery &&
                    find(vecNodePeriphery.begin(), vecNodePeriphery.end(), *it_i) == vecNodePeriphery.end())
                    vecNodePeriphery.push_back(*it_i);
            }
        }

        // cout << "*CONFIRMED* number of candidate nodes at colony periphery : " << vecNodePeriphery.size() << endl;
    }
}
void sortClonePeriphery(vector<Node> &vecNode, vector<int> &vecNodePeriphery)
{
    // function: sort according to angular locations
    unordered_map<double, int> mapAngleNode;
    vector<double> vecAngle;
    for (vector<int>::iterator it_i = vecNodePeriphery.begin(); it_i != vecNodePeriphery.end(); it_i++)
    {
        dVec coord = vecNode[*it_i].get_coord();
        double ang = atan2(coord.y, coord.x);
        while (find(vecAngle.begin(), vecAngle.end(), ang) != vecAngle.end())
            ang += 0.01 * PI;
        vecAngle.push_back(ang);
        mapAngleNode[ang] = *it_i;
    }

    // sort the vector of angles in ascending order
    sort(vecAngle.begin(), vecAngle.end());

    // clear vecNodePeriphery as information is transfered to mapAngleNode
    vecNodePeriphery.clear();
    for (vector<double>::iterator it_d = vecAngle.begin(); it_d != vecAngle.end(); it_d++)
    {
        vecNodePeriphery.push_back(mapAngleNode[*it_d]);
    }
}
void fillGapPeriphery(vector<Node> &vecNode, vector<int> &vecNodePeriphery)
{
    // function:
    // ... by checking (1) if two consecutive peripheral cells are NOT in contact and (2) if they have common neighbor cells
    // ... add their common neighbor to the vector

    vector<int> vecNodePeriphery_temp;

    for (vector<int>::iterator it_i = vecNodePeriphery.begin(); it_i != vecNodePeriphery.end(); it_i++)
    {
        int nId1 = *it_i, nId2;
        vecNodePeriphery_temp.push_back(nId1);

        // node2
        if (it_i - vecNodePeriphery.end() == -1)
            nId2 = vecNodePeriphery[0];
        else
            nId2 = *(it_i + 1);

        // find if nId2 is in nId1 neighbor list
        vector<int> vecNbNodeIds1 = vecNode[nId1].get_vecNbNodeIds();
        if (find(vecNbNodeIds1.begin(), vecNbNodeIds1.end(), nId2) != vecNbNodeIds1.end()) // found !
            continue;
        else // not found
        {
            // try to find if they share common neighbors
            vector<int> vecNbNodeIds2 = vecNode[nId2].get_vecNbNodeIds();
            for (vector<int>::iterator it_j = vecNbNodeIds1.begin(); it_j != vecNbNodeIds1.end(); it_j++)
            {
                if (find(vecNbNodeIds2.begin(), vecNbNodeIds2.end(), *it_j) != vecNbNodeIds2.end() &&
                    find(vecNodePeriphery.begin(), vecNodePeriphery.end(), *it_j) == vecNodePeriphery.end() &&
                    find(vecNodePeriphery_temp.begin(), vecNodePeriphery_temp.end(), *it_j) == vecNodePeriphery_temp.end())
                {
                    vecNodePeriphery_temp.push_back(*it_j);
                    break; // only include 1 common neighbor
                }
            }
        }
    }

    // cout << ">>Number of peripheral cells<<" << endl;
    // cout << "before filling gap: " << vecNodePeriphery.size() << endl;
    // cout << "after  filling gap: " << vecNodePeriphery_temp.size() << endl;
    // cout << "   net add:         " << vecNodePeriphery_temp.size() - vecNodePeriphery.size() << endl;

    vecNodePeriphery.clear();
    vecNodePeriphery = vecNodePeriphery_temp;
}
void calcCurvature(vector<Node> &vecNode, vector<int> &vecNodePeriphery, unordered_map<int, double> &mapNodeCurvature)
{

    double x1, y1, x2, y2, x3, y3;
    double A, B, C, D, r, kappa;
    bool flagIsConvex = true;
    for (vector<int>::iterator it_i = vecNodePeriphery.begin(); it_i != vecNodePeriphery.end(); it_i++)
    {
        x2 = vecNode[*it_i].get_coord().x;
        y2 = vecNode[*it_i].get_coord().y;

        // x1, y1
        if (it_i - vecNodePeriphery.begin() == 0)
        {
            x1 = vecNode[vecNodePeriphery[vecNodePeriphery.size() - 1]].get_coord().x;
            y1 = vecNode[vecNodePeriphery[vecNodePeriphery.size() - 1]].get_coord().y;
        }
        else
        {
            x1 = vecNode[*(it_i - 1)].get_coord().x;
            y1 = vecNode[*(it_i - 1)].get_coord().y;
        }
        // x3, y3
        if (it_i - vecNodePeriphery.end() == -1)
        {
            x3 = vecNode[vecNodePeriphery[0]].get_coord().x;
            y3 = vecNode[vecNodePeriphery[0]].get_coord().y;
        }
        else
        {
            x3 = vecNode[*(it_i + 1)].get_coord().x;
            y3 = vecNode[*(it_i + 1)].get_coord().y;
        }

        // concavity, convexity
        dPair vec21 = {x1 - x2, y1 - y2};
        dPair vec23 = {x3 - x2, y3 - y2};
        double ang_21_23 = atan2(vec21.x * vec23.y - vec21.y * vec23.x, vec21.x * vec23.x + vec21.y * vec23.y);
        if (ang_21_23 <= 0)
            flagIsConvex = true;
        else
            flagIsConvex = false;

        // A, B, C, D
        A = x1 * (y2 - y3) - y1 * (x2 - x3) + x2 * y3 - x3 * y2;
        B = (x1 * x1 + y1 * y1) * (y3 - y2) + (x2 * x2 + y2 * y2) * (y1 - y3) + (x3 * x3 + y3 * y3) * (y2 - y1);
        C = (x1 * x1 + y1 * y1) * (x2 - x3) + (x2 * x2 + y2 * y2) * (x3 - x1) + (x3 * x3 + y3 * y3) * (x1 - x2);
        D = (x1 * x1 + y1 * y1) * (x3 * y2 - x2 * y3) + (x2 * x2 + y2 * y2) * (x1 * y3 - x3 * y1) + (x3 * x3 + y3 * y3) * (x2 * y1 - x1 * y2);

        // r, kappa
        r = sqrt((B * B + C * C - 4 * A * D) / (4 * A * A));
        kappa = 1. / r;

        if (flagIsConvex == true)
        {
            // cout << "nId = " << *it_i << "; kappa = " << kappa << endl;
            mapNodeCurvature[*it_i] = kappa;
        }
        if (flagIsConvex == false)
        {
            // cout << "nId = " << *it_i << "; kappa = " << -kappa << endl;
            mapNodeCurvature[*it_i] = -kappa;
        }
    }
}

// function printGroupInfo() is an iterative function
double printGroupInfo(vector<Group> &vecGroup, vector<Node> &vecNode, double t_now, int seed)
{
    double x_min = 0, x_max = 0, y_min = 0, y_max = 0, z_min = 0, z_max = 0;
    int colony_size = 0;
    for (vector<Group>::iterator it_c = vecGroup.begin(); it_c != vecGroup.end(); it_c++)
    {
        vector<int> vecBeadIds = it_c->get_vecNodeIds();
        int groupId = it_c->get_id();
        cout << "Group_id: " << groupId << "; Group_size: " << vecBeadIds.size() << " cells" << endl;
        colony_size += vecBeadIds.size();

        int beadIdM;
        dVec coordM, velocM;
        double speedSqM = -1.;
        double speedSqAvg = 0;
        double speedAvg = 0;

        double pressureAvg = 0;
        double pressureMax = -1000;
        double pressureMin = 1000;
        int beadIdMaxP, beadIdMinP;

        for (vector<int>::iterator it_n = vecBeadIds.begin(); it_n != vecBeadIds.end(); it_n++)
        {
            int beadId = *it_n;
            // check if the cell isAlive
            bool isAlive = vecNode[beadId].get_isAlive();
            if (isAlive)
            {
                dVec veloc = vecNode[beadId].get_veloc();
                dVec coord = vecNode[beadId].get_coord();
                double speedSq = veloc.x * veloc.x + veloc.y * veloc.y + veloc.z * veloc.z;
                // speedSqAvg += speedSq;
                speedAvg += sqrt(speedSq);
                if (speedSq > speedSqM)
                {
                    speedSqM = speedSq;
                    beadIdM = beadId;
                    coordM = coord;
                    velocM = veloc;
                }

                double pressure = vecNode[beadId].get_pressure();
                pressureAvg += pressure;
                if (pressure < pressureMin)
                {
                    pressureMin = pressure;
                    beadIdMinP = beadId;
                }
                if (pressure > pressureMax)
                {
                    pressureMax = pressure;
                    beadIdMaxP = beadId;
                }

                // if (pressure < 0)
                //     cout << beadId << " pressure : " << pressure << endl;

                if (coord.x < x_min)
                    x_min = coord.x;
                if (coord.x > x_max)
                    x_max = coord.x;
                if (coord.y < y_min)
                    y_min = coord.y;
                if (coord.y > y_max)
                    y_max = coord.y;
                if (coord.z < z_min)
                    z_min = coord.z;
                if (coord.z > z_max)
                    z_max = coord.z;
            }
        }
    }

    cout << "--- Summary of Virtual Tumour ---" << endl;
    cout << "DIM_X (um) \t DIM_Y (um) \t DIM_Z (um) \t NUM_CELL" << endl;
    cout << x_max - x_min << "\t" << y_max - y_min << "\t" << z_max - z_min << "\t" << colony_size << endl;

    // write tumour size into file
    ofstream tumourFile;
    string tumourFileName;
    stringstream ssTime;
    ssTime << fixed << setprecision(1) << t_now;

    stringstream PROC_ID_SS;
    PROC_ID_SS << PROC_ID;

    stringstream seed_SS;
    seed_SS << seed;

    tumourFileName = "PID_" + PROC_ID_SS.str() + "_seed_" + seed_SS.str() + "_tumour_size.txt";
    tumourFile.open(tumourFileName, ios::app | ios::binary);
    if (t_now == 0)
        tumourFile << "time\txdim\tydim\tncell" << endl; // NEED TO WRITE NUMBER OF APOPTOTIC CELLS!!!!
    stringstream ssxdim, ssydim, ssnc;
    ssxdim << fixed << setprecision(3) << x_max - x_min;
    ssydim << fixed << setprecision(3) << y_max - y_min;
    ssnc << fixed << setprecision(1) << vecNode.size();
    tumourFile << ssTime.str() << "\t" << ssxdim.str() << "\t" << ssydim.str() << "\t" << ssnc.str() << endl;
    tumourFile.close();

    return 0.5 * (x_max - x_min + y_max - y_min);
}

// [older way of saving outputs in PDB format for visualisation using PyMol] function writeGroupInfo() is an iterative function
void writeGroupInfo(vector<Group> &vecGroup, vector<Node> &vecNode, double t_now, vector<Node> &vecNodeBd)
{
    // write in PDB format
    // ....... format ........ see https://www.ichemlabs.com/166
    // ATOM   {atom_id as global bead id} CA  ASP {chain_id as A, B, C, ...}  {residue_id as bead_id in Chromosome} {x} {y} {z} 1.00 0.00
    // TER  {atom_id_prev + 1}  {empty} ASP {chain_id}  {residue_id}
    bool flagWritePDB = true;
    if (flagWritePDB)
    {
        // get clone center of mass
        vector<int> vecNodePeriphery;
        dVec cloneCOM = {0, 0, 0};
        int nNode = vecNode.size();
        bool flagCheckPeriphery = true;
        currCloneCenterPeriphery(vecNode, vecNodePeriphery, cloneCOM, flagCheckPeriphery);
        // cout << cloneCOM.x << ", " << cloneCOM.y << ", " << cloneCOM.z << endl;
        unordered_map<int, double> mapNodeCurvature;
        sortClonePeriphery(vecNode, vecNodePeriphery);
        fillGapPeriphery(vecNode, vecNodePeriphery);
        calcCurvature(vecNode, vecNodePeriphery, mapNodeCurvature);
        // file to write coordinates of peripheral cells
        ofstream periFile;
        string periFileName;
        stringstream PROC_ID_SS;
        PROC_ID_SS << PROC_ID;
        periFileName = "PID_" + PROC_ID_SS.str() + "_ColonyPeriphery.txt";
        periFile.open(periFileName, ios::app | ios::binary);
        if (t_now == 0)
            periFile << "time\tColonySize\tCellID\tx\ty" << endl;
        for (vector<int>::iterator it_i = vecNodePeriphery.begin(); it_i != vecNodePeriphery.end(); it_i++)
        {
            stringstream sst, sscs, ssid, ssx, ssy;
            sst << t_now;
            sscs << vecNode.size();
            ssid << *it_i;
            ssx << vecNode[*it_i].get_coord().x;
            ssy << vecNode[*it_i].get_coord().y;
            // ssx << vecNode[*it_i].get_coord_prev().x;
            // ssy << vecNode[*it_i].get_coord_prev().y;
            periFile << sst.str() << "\t" << sscs.str() << "\t" << ssid.str() << "\t"
                     << ssx.str() << "\t" << ssy.str() << endl;
        }
        periFile.close();

        // file to write cell color information
        // header: node_id, group_chain, node_id_in_group, cell_generation
        ofstream nodeFile;
        string nodeFileName;
        stringstream ssPDB;
        ssPDB << fixed << setprecision(1) << t_now;
        nodeFileName = "PID_" + PROC_ID_SS.str() + "_nodeInfoColor_" + ssPDB.str() + ".txt";
        nodeFile.open(nodeFileName, ios::app | ios::binary);
        nodeFile << "CellID\tGroupLabel\tCellIDinGroup\tGeneration\tPressure(Pa)\tDensity1\tDensity2\tDensity3\tx(um)\ty(um)\tForceRad(pN)\tForceAzi(pN)\tProlifOn\tTdIntrinsic\tIsPeriphery\tnumNeighbor\tmaxCompression\tcurvature" << endl;

        // pdb file
        ofstream nodeFilePDB;
        string nodeFileNamePDB;
        // stringstream ssPDB;
        // ssPDB << fixed << setprecision(1) << t_now;
        // stringstream PROC_ID_SS; PROC_ID_SS << PROC_ID;
        nodeFileNamePDB = "PID_" + PROC_ID_SS.str() + "_tumourPDB_" + ssPDB.str() + ".pdb";

        nodeFilePDB.open(nodeFileNamePDB, ios::app | ios::binary);
        int nGroup = vecGroup.size();
        unordered_map<int, string> int2letter;
        int2letter[1] = "A";
        int2letter[2] = "B";
        int2letter[3] = "C";
        int2letter[4] = "D";
        int2letter[5] = "E";
        int2letter[6] = "F";
        int2letter[7] = "G";
        int2letter[8] = "H";
        int2letter[9] = "I";
        int2letter[10] = "J";
        int2letter[11] = "K";
        int2letter[12] = "L";
        int2letter[13] = "M";
        int2letter[14] = "N";
        int2letter[15] = "O";
        int2letter[16] = "P";
        int2letter[17] = "Q";
        int2letter[18] = "R";
        int2letter[19] = "S";
        int2letter[20] = "T";
        int2letter[21] = "U";
        int2letter[22] = "V";
        int2letter[23] = "W";
        int2letter[24] = "X";
        int2letter[25] = "Y";
        int2letter[26] = "Z";
        // add more groups to "Z"
        for (int more = 27; more <= 200; more++)
            int2letter[more] = "Z";

        // version 1: write group by group
        // ... problem: the global atom ids don't follow sequential order
        if (false)
        {
            int cntAll = 0;
            for (vector<Group>::iterator it_c = vecGroup.begin(); it_c != vecGroup.end(); it_c++)
            {
                vector<int> vecNodeIds = it_c->get_vecNodeIds();
                int groupId = it_c->get_id();

                int groupSegSize = 9999; // this is used for visualization purpose when the length is too long
                int segId = 0;
                for (vector<int>::iterator it_n = vecNodeIds.begin(); it_n != vecNodeIds.end(); it_n++)
                {
                    int nodeId = *it_n;
                    dVec coord = vecNode[nodeId].get_coord();
                    // dVec coord = vecNode[nodeId].get_coord_prev();
                    double x = coord.x, y = coord.y, z = coord.z;

                    stringstream ss0, ss1;
                    stringstream ssx, ssy, ssz;
                    ss0 << setw(5) << nodeId + groupId;
                    ss1 << setw(4) << 1 + (it_n - vecNodeIds.begin());
                    if (x < -1000)
                        ssx << fixed << setw(8) << setprecision(2) << x;
                    else
                        ssx << fixed << setw(8) << setprecision(3) << x;
                    if (y < -1000)
                        ssy << fixed << setw(8) << setprecision(2) << y;
                    else
                        ssy << fixed << setw(8) << setprecision(3) << y;
                    if (z < -1000)
                        ssz << fixed << setw(8) << setprecision(2) << z;
                    else
                        ssz << fixed << setw(8) << setprecision(3) << z;

                    int nodeRank = (it_n - vecNodeIds.begin());
                    segId = nodeRank / groupSegSize;
                    stringstream ss0_new, ss1_new;
                    ss0_new << setw(5) << nodeId + groupId + segId;
                    ss1_new << setw(4) << 1 + nodeRank % groupSegSize;

                    if (nodeRank >= groupSegSize && nodeRank % groupSegSize == 0) // break the Group into "segments"
                    {
                        stringstream ss2;
                        ss2 << setw(5) << nodeRank + groupId;
                        nodeFilePDB << "TER   " << ss2.str() << "    " << "  " << "ASP" << " " << int2letter[groupId] << endl;
                    }

                    // the following line uses groupId+segId before
                    nodeFilePDB << "ATOM  " << ss0_new.str() << "  CA" << "  " << "ASP" << " " << int2letter[groupId] << ss1_new.str() << "    "
                                << ssx.str() << ssy.str() << ssz.str() << "  1.00" << "  0.00" << endl;

                    // write into file containing color information
                    int gen = vecNode[nodeId].get_gen();
                    double pressure = vecNode[nodeId].get_pressure();
                    unordered_map<string, int> umapCD = vecNode[nodeId].get_umapCellDensity();
                    double time_double_intrinsic = vecNode[nodeId].get_time_double_intrinsic();
                    double time_double = vecNode[nodeId].get_time_double();
                    short flagProlifOn = 1, flagIsPeriphery = 0;
                    if (time_double < 0)
                        flagProlifOn = 0;
                    if (find(vecNodePeriphery.begin(), vecNodePeriphery.end(), nodeId) != vecNodePeriphery.end())
                        flagIsPeriphery = 1;
                    vector<int> vecNbNodeIds = vecNode[nodeId].get_vecNbNodeIds();
                    int numNeighbor = vecNbNodeIds.size();

                    dVec veloc = vecNode[nodeId].get_veloc(), unitVelocVec, zeroVec = {0, 0, 0};
                    dVec vecPos = {coord.x - cloneCOM.x, coord.y - cloneCOM.y, coord.z - cloneCOM.z};
                    double veloc_size = dVecDist(zeroVec, veloc);
                    // unitVelocVec = {veloc.x/veloc_size, veloc.y/veloc_size, veloc.z/veloc_size};
                    // double d2cloneCOM = dVecDist(cloneCOM, coord);
                    // dVec unitVec = { (coord.x-cloneCOM.x)/d2cloneCOM,
                    //                  (coord.y-cloneCOM.y)/d2cloneCOM,
                    //                  (coord.z-cloneCOM.z)/d2cloneCOM };
                    // double dot = unitVec.x*unitVelocVec.y + unitVec.y*unitVelocVec.x;
                    // double det = unitVec.x*unitVelocVec.y - unitVec.y*unitVelocVec.x;
                    double dot = vecPos.x * veloc.x + vecPos.y * veloc.y;
                    double det = vecPos.x * veloc.y - vecPos.y * veloc.x;
                    double angle = atan2(det, dot);
                    double veloc_rad = veloc_size * cos(angle), veloc_azi = veloc_size * sin(angle);
                    double factor = 1 * GAMMA_VISCOUS / 3.6 / 3.6; // speed to pN
                    double force_rad = veloc_rad * factor, force_azi = veloc_azi * factor;

                    stringstream ssg, ssp, sscd1, sscd2, sscd3, ssfr, ssfa, sspo, sstd, ssip, ssnnb;
                    ssg << gen;
                    ssp << pressure;
                    sscd1 << umapCD["d1"];
                    sscd2 << umapCD["d2"];
                    sscd3 << umapCD["d3"];
                    ssfr << force_rad;
                    ssfa << force_azi;
                    sspo << flagProlifOn;
                    sstd << time_double_intrinsic;
                    ssip << flagIsPeriphery;
                    ssnnb << numNeighbor;
                    // cout << umapCD["d1"] << ", " << umapCD["d2"] << ", " << umapCD["d3"] << endl;
                    nodeFile << ss0_new.str() << "\t" << int2letter[groupId] << "\t"
                             << ss1_new.str() << "\t" << ssg.str() << "\t" << ssp.str() << "\t"
                             << sscd1.str() << "\t" << sscd2.str() << "\t" << sscd3.str() << "\t"
                             << ssx.str() << "\t" << ssy.str() << "\t"
                             << ssfr.str() << "\t" << ssfa.str() << "\t"
                             << sspo.str() << "\t" << sstd.str() << "\t"
                             << ssip.str() << "\t" << ssnnb.str() << endl;
                }

                stringstream ss2;
                cntAll = vecNodeIds[vecNodeIds.size() - 1] + 1 + groupId;
                ss2 << setw(5) << vecNodeIds[vecNodeIds.size() - 1] + 1 + groupId;
                nodeFilePDB << "TER   " << ss2.str() << "    " << "  " << "ASP" << " " << int2letter[groupId] << endl;
            }
        }

        // version 2: write node by node
        // ... neglect: the rank of nodes in group
        if (true)
        {
            for (vector<Node>::iterator it_n = vecNode.begin(); it_n != vecNode.end(); it_n++)
            {
                // check if the cell isAlive
                bool isAlive = it_n->get_isAlive();
                int nodeId = it_n->get_id();
                int groupId = it_n->get_groupId();
                if (isAlive == false && groupId != 26)
                {
                    cout << "WARNING! this cell isAlive = " << isAlive << "; groupId = " << groupId << endl;
                }

                dVec coord = it_n->get_coord();
                double x = coord.x, y = coord.y, z = coord.z;

                stringstream ss0, ss1;
                stringstream ssx, ssy, ssz;
                ss0 << setw(7) << nodeId + 1;
                ss1 << setw(4) << 1;
                if (x < -1000)
                    ssx << fixed << setw(8) << setprecision(2) << x;
                else
                    ssx << fixed << setw(8) << setprecision(3) << x;
                if (y < -1000)
                    ssy << fixed << setw(8) << setprecision(2) << y;
                else
                    ssy << fixed << setw(8) << setprecision(3) << y;
                if (z < -1000)
                    ssz << fixed << setw(8) << setprecision(2) << z;
                else
                    ssz << fixed << setw(8) << setprecision(3) << z;

                nodeFilePDB << "ATOM" << ss0.str() << "  CA" << "  " << "ASP" << " " << int2letter[groupId] << ss1.str() << "    "
                            << ssx.str() << ssy.str() << ssz.str() << "  1.00" << "  0.00" << endl;

                // write into file containing color information
                int gen = vecNode[nodeId].get_gen();
                double pressure = vecNode[nodeId].get_pressure();
                unordered_map<string, int> umapCD = vecNode[nodeId].get_umapCellDensity();
                double time_double_intrinsic = vecNode[nodeId].get_time_double_intrinsic();
                double time_double = vecNode[nodeId].get_time_double();
                short flagProlifOn = 1, flagIsPeriphery = 0;
                double curvature = 0;
                if (time_double < 0)
                    flagProlifOn = 0;
                if (find(vecNodePeriphery.begin(), vecNodePeriphery.end(), nodeId) != vecNodePeriphery.end())
                {
                    flagIsPeriphery = 1;
                    if (vecNodePeriphery.size() > 10)
                        curvature = mapNodeCurvature[nodeId];
                }
                vector<int> vecNbNodeIds = vecNode[nodeId].get_vecNbNodeIds();
                int numNeighbor = vecNbNodeIds.size();

                dVec veloc = vecNode[nodeId].get_veloc(), unitVelocVec, zeroVec = {0, 0, 0};
                dVec vecPos = {coord.x - cloneCOM.x, coord.y - cloneCOM.y, coord.z - cloneCOM.z};
                double veloc_size = dVecDist(zeroVec, veloc);

                double dot = vecPos.x * veloc.x + vecPos.y * veloc.y;
                double det = vecPos.x * veloc.y - vecPos.y * veloc.x;
                double angle = atan2(det, dot);
                double veloc_rad = veloc_size * cos(angle), veloc_azi = veloc_size * sin(angle);
                double factor = 1 * GAMMA_VISCOUS / 3.6 / 3.6; // speed to pN
                double force_rad = veloc_rad * factor, force_azi = veloc_azi * factor;

                vector<double> vecNbPressures = vecNode[nodeId].get_vecNbPressures();
                double maxCompr = 0;
                if (vecNbPressures.size())
                    maxCompr = *min_element(vecNbPressures.begin(), vecNbPressures.end());

                stringstream ssg, ssp, sscd1, sscd2, sscd3, ssfr, ssfa, sspo, sstd, ssip, ssnnb, ssmc, sscur;
                ssg << gen;
                ssp << pressure;
                sscd1 << umapCD["d1"];
                sscd2 << umapCD["d2"];
                sscd3 << umapCD["d3"];
                ssfr << force_rad;
                ssfa << force_azi;
                sspo << flagProlifOn;
                sstd << time_double_intrinsic;
                ssip << flagIsPeriphery;
                ssnnb << numNeighbor;
                ssmc << maxCompr;
                sscur << curvature;

                // cout << umapCD["d1"] << ", " << umapCD["d2"] << ", " << umapCD["d3"] << endl;
                nodeFile << ss0.str() << "\t" << int2letter[groupId] << "\t"
                         << ss1.str() << "\t" << ssg.str() << "\t" << ssp.str() << "\t"
                         << sscd1.str() << "\t" << sscd2.str() << "\t" << sscd3.str() << "\t"
                         << ssx.str() << "\t" << ssy.str() << "\t"
                         << ssfr.str() << "\t" << ssfa.str() << "\t"
                         << sspo.str() << "\t" << sstd.str() << "\t"
                         << ssip.str() << "\t" << ssnnb.str() << "\t"
                         << ssmc.str() << "\t" << sscur.str() << endl;
            }

            // update: also write vecNodeBd
            const bool flagWriteNodeBd = true;
            if (flagWriteNodeBd)
            {
                for (vector<Node>::iterator it_n_bd = vecNodeBd.begin(); it_n_bd != vecNodeBd.end(); it_n_bd++)
                {
                    dVec coord_bd = it_n_bd->get_coord();
                    double x_bd = coord_bd.x, y_bd = coord_bd.y, z_bd = coord_bd.z;
                    int nodeId_bd = vecNode.size() + it_n_bd->get_id();

                    stringstream ss0_bd, ss1_bd;
                    stringstream ssx_bd, ssy_bd, ssz_bd;
                    ss0_bd << setw(7) << nodeId_bd + 1;
                    ss1_bd << setw(4) << 1;
                    if (x_bd < -1000)
                        ssx_bd << fixed << setw(8) << setprecision(2) << x_bd;
                    else
                        ssx_bd << fixed << setw(8) << setprecision(3) << x_bd;
                    if (y_bd < -1000)
                        ssy_bd << fixed << setw(8) << setprecision(2) << y_bd;
                    else
                        ssy_bd << fixed << setw(8) << setprecision(3) << y_bd;
                    if (z_bd < -1000)
                        ssz_bd << fixed << setw(8) << setprecision(2) << z_bd;
                    else
                        ssz_bd << fixed << setw(8) << setprecision(3) << z_bd;

                    nodeFilePDB << "ATOM" << ss0_bd.str() << "  CA" << "  " << "ASP" << " " << "W" << ss1_bd.str() << "    "
                                << ssx_bd.str() << ssy_bd.str() << ssz_bd.str() << "  1.00" << "  0.00" << endl;
                }
            }

            stringstream ss2;
            int cntAll = vecNode.size() + 1;
            if (flagWriteNodeBd)
                cntAll += vecNodeBd.size();
            ss2 << setw(5) << cntAll;
            nodeFilePDB << "TER   " << ss2.str() << "    " << "  " << "ASP" << " " << int2letter[vecGroup.size() + 1] << endl;
        }

        nodeFile.close();
        nodeFilePDB.close();
    }
}

// function writeNodeDynamics() is an iterative function
void writeNodeDynamics(vector<Node> &vecNode, double t_now, int seed, bool flagWriteMut)
{
    ofstream nodeFile;
    string nodeFileName;
    stringstream PROC_ID_SS;
    PROC_ID_SS << PROC_ID;

    stringstream seed_SS;
    seed_SS << seed;

    nodeFileName = "PID_" + PROC_ID_SS.str() + "_seed_" + seed_SS.str() + "_nodeDynamics.txt";

    nodeFile.open(nodeFileName, ios::app | ios::binary);
    if (t_now == 0)
        nodeFile << "t(h)\tCellID\tLineageID\tGroupID\tgeneration\tx(um)\ty\tz" << endl;

    // [NOT USED IN THE LTM STUDY] mutation file 1 : record the gene mutations in each cell
    // header: t, id, mutation_ids
    // ofstream nodeMutFile;
    // string nodeMutFileName;
    // nodeMutFileName = "PID_" + PROC_ID_SS.str() + "_nodeMutations.txt";
    // nodeMutFile.open(nodeMutFileName, ios::app | ios::binary);
    // if (t_now == 0)
    //     nodeMutFile << "t(h)\tCellID\tMutations" << endl;

    // [NOT USED IN THE LTM STUDY] mutation file 2 : record cell ids that harbor driver genes
    // header: t, driver_id, cell_ids
    // ofstream nodeMutFile2;
    // string nodeMutFileName2;
    // nodeMutFileName2 = "PID_" + PROC_ID_SS.str() + "_nodeHarborDrivers.txt";
    // nodeMutFile2.open(nodeMutFileName2, ios::app | ios::binary);
    // if (t_now == 0)
    //     // nodeFile << "t(h)\tCellID\tLineageID\tGroupID\tgeneration\tx(um)\ty\tz\tvx(um/s)\tvy\tvz" << endl;
    //     nodeMutFile2 << "t(h)\tMutationID\tCellIds" << endl;
    // unordered_map<int, vector<int>> umapDriverCells; // driver mutation id : vector of cell ids
    // for (int driver = 0; driver < NUM_MUT_DRIVER; driver++)
    // {
    //     vector<int> vecCellsHarborDriver;
    //     umapDriverCells[driver] = vecCellsHarborDriver;
    // }

    // [NOT USED IN THE LTM STUDY] mutation file 3 : record number of cells harboring a given gene
    // header: t, mutation_id, number_of_cells
    // ofstream nodeMutFile3;
    // string nodeMutFileName3;
    // nodeMutFileName3 = "PID_" + PROC_ID_SS.str() + "_mutationCellNumber.txt";
    // nodeMutFile3.open(nodeMutFileName3, ios::app | ios::binary);
    // if (t_now == 0)
    //     nodeMutFile3 << "t(h)\tMutationID\tCellNumber" << endl;
    // unordered_map<int, int> umapMutationCellNumber; // mutation id : cell number
    // for (int mut = 0; mut < NUM_MUT_POOL; mut++)
    //     umapMutationCellNumber[mut] = 0;

    for (vector<Node>::iterator it_s = vecNode.begin(); it_s != vecNode.end(); it_s++)
    {
        // check if the cell isAlive
        bool isAlive = it_s->get_isAlive();
        if (isAlive)
        {
            // write node dynamics
            int nId = it_s->get_id();
            int lId = it_s->get_lineageId();
            int gId = it_s->get_groupId();
            int gen = it_s->get_gen();
            dVec coord = it_s->get_coord(), veloc = it_s->get_veloc();
            stringstream sst, ssid, sslid, ssgid, ssg, ssx, ssy, ssz, ssvx, ssvy, ssvz;
            sst << t_now;
            ssid << nId;
            sslid << lId;
            ssgid << gId;
            ssg << gen;
            ssx << coord.x;
            ssy << coord.y;
            ssz << coord.z;
            ssvx << veloc.x;
            ssvy << veloc.y;
            ssvz << veloc.z;

            nodeFile << sst.str() << "\t" << ssid.str() << "\t" << sslid.str() << "\t" << ssgid.str() << "\t"
                     << ssg.str() << "\t"
                     << ssx.str() << "\t" << ssy.str() << "\t" << ssz.str() << "\t" << endl;

            // [NOT USED IN THE LTM STUDY] write mutation ids
            // if (TYPE_SUBCL == "emerge")
            // {
            //     if (flagWriteMut)
            //         nodeMutFile << sst.str() << "\t" << ssid.str() << "\t";
            //     vector<int> vecMutIds = it_s->get_vecMutIds();
            //     for (vector<int>::iterator it_m = vecMutIds.begin(); it_m != vecMutIds.end(); it_m++)
            //     {
            //         if (flagWriteMut)
            //         {
            //             stringstream ssm;
            //             ssm << *(it_m);
            //             nodeMutFile << ssm.str() << "\t";
            //         }

            //         // if this is a driver mutation, attach the cell to the umap
            //         if (*(it_m) < NUM_MUT_DRIVER)
            //             umapDriverCells[*(it_m)].push_back(nId);

            //         // add the cell count to umapMutationCellNumber
            //         umapMutationCellNumber[*(it_m)] += 1;
            //     }
            //     if (flagWriteMut)
            //         nodeMutFile << endl;
            // }
        }
    }
    nodeFile.close();
    // nodeMutFile.close();

    // [NOT USED IN THE LTM STUDY]
    // if (TYPE_SUBCL == "emerge")
    // {
    //     // write to mutation file 2
    //     for (unordered_map<int, vector<int>>::iterator it_um = umapDriverCells.begin(); it_um != umapDriverCells.end(); it_um++)
    //     {
    //         stringstream sst, ssdriver;
    //         int driver = it_um->first;
    //         vector<int> vecCellsHarborDriver = it_um->second;

    //         sst << t_now;
    //         ssdriver << driver;
    //         nodeMutFile2 << sst.str() << "\t" << ssdriver.str() << "\t";

    //         for (vector<int>::iterator it_n = vecCellsHarborDriver.begin(); it_n != vecCellsHarborDriver.end(); it_n++)
    //         {
    //             stringstream ssid;
    //             ssid << *(it_n);
    //             nodeMutFile2 << ssid.str() << "\t";
    //         }
    //         nodeMutFile2 << endl;
    //     }

    //     // write to mutation file 3
    //     for (unordered_map<int, int>::iterator it_um = umapMutationCellNumber.begin(); it_um != umapMutationCellNumber.end(); it_um++)
    //     {
    //         stringstream sst, ssmut, sscnum;
    //         int mut = it_um->first;
    //         int cnum = it_um->second;

    //         sst << t_now;
    //         ssmut << mut;
    //         sscnum << cnum;

    //         nodeMutFile3 << sst.str() << "\t" << ssmut.str() << "\t" << sscnum.str() << endl;
    //     }
    // }
    // nodeMutFile2.close();
    // nodeMutFile3.close();
}
// [NOT USED IN THE LTM STUDY] function writeNodeMitosis() is an iterative function
// ... but the information is only recorded if a cell isAlive and undergoes mitosis
// ... so no need to check if a cell isAlive
void writeNodeMitosis(vector<unordered_map<string, double>> &mitosisToWrite, double t_now)
{
    // header: cloneSize, cloneRad, nodeId, d2CloneCenter, d2ClonePeriphery
    vector<string> vecInfo;
    vecInfo.push_back("cloneSize");
    vecInfo.push_back("cloneRad");
    vecInfo.push_back("nodeId");
    vecInfo.push_back("d2CloneCenter");
    vecInfo.push_back("d2ClonePeriphery");
    vecInfo.push_back("time_double_intrinsic");
    vecInfo.push_back("time_double_actual");
    vecInfo.push_back("generation");
    vecInfo.push_back("nodeId_child");
    vecInfo.push_back("lineageId");

    ofstream nodeFile;
    string nodeFileName;
    stringstream ss;
    ss << fixed << setprecision(1) << t_now;
    stringstream PROC_ID_SS;
    PROC_ID_SS << PROC_ID;
    nodeFileName = "PID_" + PROC_ID_SS.str() + "_nodeMitosis_" + ss.str() + ".txt";
    nodeFile.open(nodeFileName, ios::app | ios::binary);

    // header
    for (vector<string>::iterator it_s = vecInfo.begin(); it_s != vecInfo.end(); it_s++)
    {
        string info = *it_s;
        nodeFile << info << "\t";
    }
    nodeFile << endl;

    // content
    for (vector<unordered_map<string, double>>::iterator it_um = mitosisToWrite.begin(); it_um != mitosisToWrite.end(); it_um++)
    {
        unordered_map<string, double> mapInfo = *it_um;
        for (vector<string>::iterator it_s = vecInfo.begin(); it_s != vecInfo.end(); it_s++)
        {
            string info = *it_s;
            stringstream ss_info;
            ss_info << fixed << setprecision(2) << mapInfo[info];
            nodeFile << ss_info.str() << "\t";
        }
        nodeFile << endl;
    }
    nodeFile.close();
}

double dVecDist(dVec &p1, dVec &p2)
{
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
}
