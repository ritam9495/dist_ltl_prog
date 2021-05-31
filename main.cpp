#include "z3++.h"
#include <string>
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <time.h>
#include <omp.h>
#include <stdlib.h>

using namespace std;
using namespace z3;

class createSMTMon {
private:
    string traceFileName;
    string EValue[10][200];
    string EType[10][200];
    string EName[10][200];
    int ETime[10][200][3];
    int numProcess, numEvent[10];

    string EValueAll[10][200];
    string ETypeAll[10][200];
    int ETimeAll[10][200][3];
    int numEventAll[10];
    int hbSet[500][2];
    int numHbSet;

    void readFile();
    void readSegTrace(int, int, int);

public:
    createSMTMon();
    createSMTMon(string);
    void makeSMTMon(int, int, int, int);
    int solveSMTMon(string, string);
};

createSMTMon::createSMTMon() {
    traceFileName = "trace2_20.txt";
    int i, j;
    for(i = 0; i < 10; i++) {
        for(j = 0; j < 200; j++) {
            EName[i][j] = "";
            EValue[i][j] = "";
            EType[i][j] = "";
            ETime[i][j][0] = 0;
            ETime[i][j][1] = 0;
            ETime[i][j][2] = 0;

            EValueAll[i][j] = "";
            ETypeAll[i][j] = "";
            ETimeAll[i][j][0] = 0;
            ETimeAll[i][j][1] = 0;
            ETimeAll[i][j][2] = 0;
        }
        numEvent[i] = 0;
        numEventAll[i] = 0;
    }
    numProcess = 0;
    for(i = 0; i < 500; i++) {
        hbSet[i][0] = 0;
        hbSet[i][1] = 0;
    }
    numHbSet = 0;
}

createSMTMon::createSMTMon(string one) {
    traceFileName = one;
    int i, j;
    for(i = 0; i < 10; i++) {
        for(j = 0; j < 200; j++) {
            EName[i][j] = "";
            EValue[i][j] = "";
            EType[i][j] = "";
            ETime[i][j][0] = 0;
            ETime[i][j][1] = 0;
            ETime[i][j][2] = 0;

            EValueAll[i][j] = "";
            ETypeAll[i][j] = "";
            ETimeAll[i][j][0] = 0;
            ETimeAll[i][j][1] = 0;
            ETimeAll[i][j][2] = 0;
        }
        numEvent[i] = 0;
        numEventAll[i] = 0;
    }
    numProcess = 0;
    for(i = 0; i < 500; i++) {
        hbSet[i][0] = 0;
        hbSet[i][1] = 0;
    }
    numHbSet = 0;
}

void createSMTMon::readFile() {
    ifstream MyFile(traceFileName);
    string line;
    int proc = 0, event = -1, num = 0;
    string word;

    while(getline(MyFile, line)) {
        if(line.find("\\Process") != -1) {
            // cout << "1" << line << "\n";
            numProcess++;
        }
        else if(line.find("Event") != -1) {
            // cout << "2" << line << "\n";
            int a1 = line.find('\"');
            int a2 = line.find('\"', a1 + 1);
            int b1 = line.find('\"', a2 + 1);
            int b2 = line.find('\"', b1 + 1);
            int c1 = line.find('\"', b2 + 1);
            int c2 = line.find('\"', c1 + 1);
            int d1 = line.find('\"', c2 + 1);
            int d2 = line.find(',', d1 + 1);
            int d3 = line.find(',', d2 + 1);
            int d4 = line.find('\"', d3 + 1);
            EType[numProcess][numEventAll[numProcess]] = line.substr(b1 + 1, b2 - b1 - 1);
            EValue[numProcess][numEventAll[numProcess]] = line.substr(c1 + 1, c2 - c1 - 1);
            ETimeAll[numProcess][numEventAll[numProcess]][0] = stoi(line.substr(d1 + 2, d2 - d1));
            ETimeAll[numProcess][numEventAll[numProcess]][1] = stoi(line.substr(d2 + 1, d3 - d2));
            ETimeAll[numProcess][numEventAll[numProcess]][2] = stoi(line.substr(d3 + 1, d4 - d3 - 1));
            // cout << EType[numProcess][numEventAll[numProcess]] << ":";
            // cout << EValue[numProcess][numEventAll[numProcess]] << ":";
            // cout << ETimeAll[numProcess][numEventAll[numProcess]][0] << ", " << ETimeAll[numProcess][numEventAll[numProcess]][1] << ", " << ETimeAll[numProcess][numEventAll[numProcess]][2] << ":";
            numEventAll[numProcess]++;
            // cout << numEventAll[numProcess] << "\n";
        }
    }
    MyFile.close();
}

void createSMTMon::readSegTrace(int start, int end, int eps) {
    int num = 1;
    start = max(start - eps, 0);
    // cout << start << " : " << end << "\n";
    for(int i = 0; i < numProcess; i++) {
        numEvent[i] = 0;
        for(int j = 0; j < numEventAll[i]; j++) {
            if(ETimeAll[i][j][0] >= start && ETimeAll[i][j][0] <= end) {
                numEvent[i]++;
                EName[i][j] = to_string(num++);
                EType[i][j] = ETypeAll[i][j];
                EValue[i][j] = EValueAll[i][j];
                ETime[i][j][0] = ETimeAll[i][j][0];
                ETime[i][j][1] = ETimeAll[i][j][1];
                ETime[i][j][2] = ETimeAll[i][j][2];
            }
        }
    }
    for(int i = 0; i < numProcess; i++) {
        for(int j = 0; j < numEvent[i]; j++) {
            if(EType[i][j] == "R") {
                for(int k = 0; k < numProcess; k ++) {
                    for(int l = 0; l < numEvent[k]; l++) {
                        if(EType[k][l] == "S" && ETime[i][j][0] == stoi(EValue[i][j].substr(0, EValue[i][j].find(",")))) {
                            hbSet[numHbSet][0] = stoi(EName[i][j]);
                            hbSet[numHbSet++][1] = stoi(EName[k][l]);
                            // cout << "Happen Before: " << hbSet[numHbSet - 1][0] << " : " << hbSet[numHbSet - 1][1] << "\n";
                        }
                    }
                }
            }
            for(int k = 0; k < numProcess; k++) {
                for(int l = 0; l < numEvent[k]; l++) {
                    if(abs(ETime[i][j][0] - ETime[k][l][0]) >= eps && i != k) {
                        hbSet[numHbSet][0] = stoi(EName[i][j]);
                        hbSet[numHbSet++][1] = stoi(EName[k][l]);
                        // cout << "Happen Before: " << hbSet[numHbSet - 1][0] << " : " << hbSet[numHbSet - 1][1] << "\n";
                    }
                }
            }
        }
    }
}

void createSMTMon::makeSMTMon(int eps, int segLength, int maxTime, int numThreads) {
    readFile();

    string statesList[10];
    string line;
    int numStates = 0;
    string currentState[10];
    int numCurrentState = 0;
    bool finalAccept = false, finalReject = false;
    // ifstream fin;
    // fin.open("formulaFiles/f" + to_str(formulaNum) + ".txt");
    // while(getline(fin, line)) {
    //     if(line.find("style=filled") != -1) {
    //         int a = (int) line.find('\"');
    //         int b = (int) line.find('\"', a + 1);
    //         statesList[numStates++] = line.substr(a + 1, b - a - 1);
    //         // cout << statesList[numStates - 1] << "\n";
    //         if(line.find("(0, 0)") != -1)
    //             currentState[numCurrentState++] = line.substr(a, b - a);
    //         else if(line.find("(-1, 1)") != -1)
    //             finalReject = true;
    //         else if(line.find("(1, -1)") != -1)
    //             finalAccept = true;
    //     }
    // }
    // fin.close();

    int resultMat[100][10][10];
    for(int  i = 0; i < int (maxTime/ segLength) + 1; i++) {
        for(int j = 0; j < numStates; j++) {
            for(int k = 0; k < numStates; k++) {
                resultMat[i][j][k] = 0;
            }
        }
    }

    int numSegment = 0;
    // numThreads = 1;
    int segmentL = maxTime/numThreads;
    // cout << "total numThreads: " << numThreads << endl;
    #pragma omp parallel num_threads(numThreads)
    {
        #pragma omp for
        for(int startT = 0; startT < maxTime; startT = startT + segmentL) {
            // cout << "Thread: " << startT << " : " << startT + segmentL << endl;

            int start, end;
            end = startT;

            while(end < startT + segmentL) {
                numSegment++;
                start = end;
                end += segLength;
                // cout << start << " : " << end << "\n";
                numStates = 1;
                readSegTrace(start, end, eps);
                for (int i = 0; i < numStates; i++) {
                    // if(statesList[i] == "(1, -1)" || statesList[i] == "(-1, 1)")
                    //     continue;
                    for (int j = 0; j < numStates; j++) {
                        // cout << statesList[i] << " : " << statesList[j] << "\n";
                        resultMat[numSegment][i][j] = solveSMTMon(statesList[i], statesList[j]);
                    }
                }
            }
        }
    }

    for(int i = 0; i < numSegment; i++) {
        string newCurrentState[10];
        int numNewCurrentState = 0;
        for(int j = 0; j < numCurrentState; j++) {
            for(int k = 0; k < numStates; k++) {
                if(statesList[k] == currentState[j]) {
                    for(int l = 0; l < numStates; l++) {
                        if(resultMat[i][k][l] == 1)
                            newCurrentState[numNewCurrentState++] = statesList[l];
                    }
                }
            }
        }
        numCurrentState = 0;
        for(int j = 0; j < numNewCurrentState; j++) {
            currentState[numCurrentState++] = newCurrentState[j];
        }
    }
}

int createSMTMon::solveSMTMon(string fromState, string toState) {
    context c;

    solver s(c);

    expr_vector eventList(c);

    int totNumEvents = 0;
    for(int i = 0; i < numProcess; i++)
        totNumEvents += numEvent[i];

    for(int i = 0; i < pow(totNumEvents, 2); i++)
    {
        std::string str = "eventList" + std::to_string(i);
        eventList.push_back(c.bv_const(str.c_str(), totNumEvents));
        s.add(eventList[i] == i);
    }

    func_decl f = z3::function("f", c.int_sort(), c.bv_sort(totNumEvents));

    s.add(f(0) == eventList[0]);

    for (int i = 0; i <= totNumEvents; i++)
    {
        expr_vector event_range(c);

        for (int j = 0; j < pow(totNumEvents, 2); j++)
        {
            event_range.push_back(f(i) == eventList[j]);
        }

        s.add(mk_or(event_range));
    }

    for(int i = 0; i < totNumEvents; i++)
    {
        expr_vector event_order(c);

        for(int j = 1; j < pow(totNumEvents, 2); j = j * 2)
            event_order.push_back(bv2int(f(i + 1), false) - bv2int(f(i), false) == j);

        // s.add(f(i + 1) > f(i));
        s.add(mk_or(event_order));
    }

    s.add(f(totNumEvents) == eventList[(int)pow(totNumEvents, 2) - 1]);

    expr x = c.int_const("x");
    s.add(0 <= x && x <= totNumEvents);

    expr_vector b1(c);
    expr_vector b2(c);

    for(int i = 0; i < numHbSet; i++)
    {
        std::string str = "b1" + std::to_string(i);
        b1.push_back(c.bv_const(str.c_str(), totNumEvents));
        str = "b2" + std::to_string(i);
        b2.push_back(c.bv_const(str.c_str(), totNumEvents));
        int p0 = (int)pow(2, hbSet[i][0] - 1);
        int p1 = (int)pow(2, hbSet[i][1] - 1);

        s.add(b1[i] == p1);
        s.add(b2[i] == p0);
        s.add(forall(x, implies(bv2int(f(x) & b1[i], false) != 0, bv2int(f(x) & b2[i], false) != 0)));
    }

    // func_decl_vector frontier(c);

    // for(int i = 0; i < numProcess; i++)
    // {
    //     std::string str = "frontier" + std::to_string(i);
    //     frontier.push_back(z3::function(str.c_str(), c.int_sort(), c.bv_sort(totNumEvents)));

    //     expr_vector frontier_order(c);

    //     for(int k = 0; k < totNumEvents; k++)
    //     {
    //         for(int j = 0; j < numEvent[i]; j++)
    //         {
    //             frontier_order.push_back(frontier[i](k) == EName[i][j]);
    //         }
    //         s.add(mk_or(frontier_order));
    //     }
    // }

    if(s.check() == sat)
    {
//        std::string verdict = std::to_string(i) + " : Sat";
        // std::cout << "Sat" << std::endl;
        return true;

//        model m = s.get_model();
//        std::cout << m << std::endl;
    }
    else
    {
//        std::string verdict = std::to_string(i) + " : Unsat";
        // std::cout << "Unsat" << std::endl;
        return false;

    }

    Z3_reset_memory();
    s.reset();
}

class createSMTProg {
private:
    string traceFileName;
    string EValue[10][200];
    string EType[10][200];
    string EName[10][200];
    int ETime[10][200][3];
    int numProcess, numEvent[10];

    string EValueAll[10][200];
    string ETypeAll[10][200];
    int ETimeAll[10][200][3];
    int numEventAll[10];
    int hbSet[500][2];
    int numHbSet;

    void readFile();
    void readSegTrace(int, int, int);

public:
    createSMTProg();
    createSMTProg(string);
    void makeSMTProg(int, int, int, string, int);
    bool solveSMTProg(string);
};

createSMTProg::createSMTProg() {
    traceFileName = "trace2_20.txt";
    int i, j;
    for(i = 0; i < 10; i++) {
        for(j = 0; j < 200; j++) {
            EName[i][j] = "";
            EValue[i][j] = "";
            EType[i][j] = "";
            ETime[i][j][0] = 0;
            ETime[i][j][1] = 0;
            ETime[i][j][2] = 0;

            EValueAll[i][j] = "";
            ETypeAll[i][j] = "";
            ETimeAll[i][j][0] = 0;
            ETimeAll[i][j][1] = 0;
            ETimeAll[i][j][2] = 0;
        }
        numEvent[i] = 0;
        numEventAll[i] = 0;
    }
    numProcess = 0;
    for(i = 0; i < 500; i++) {
        hbSet[i][0] = 0;
        hbSet[i][1] = 0;
    }
    numHbSet = 0;
}

createSMTProg::createSMTProg(string one) {
    traceFileName = one;
    int i, j;
    for(i = 0; i < 10; i++) {
        for(j = 0; j < 200; j++) {
            EName[i][j] = "";
            EValue[i][j] = "";
            EType[i][j] = "";
            ETime[i][j][0] = 0;
            ETime[i][j][1] = 0;
            ETime[i][j][2] = 0;

            EValueAll[i][j] = "";
            ETypeAll[i][j] = "";
            ETimeAll[i][j][0] = 0;
            ETimeAll[i][j][1] = 0;
            ETimeAll[i][j][2] = 0;
        }
        numEvent[i] = 0;
        numEventAll[i] = 0;
    }
    numProcess = 0;
    for(i = 0; i < 500; i++) {
        hbSet[i][0] = 0;
        hbSet[i][1] = 0;
    }
    numHbSet = 0;
}

void createSMTProg::readFile(int start, int end) {
    ifstream MyFile(traceFileName);
    string line;
    int proc = 0, event = -1, num = 0;
    string word;

    while(getline(MyFile, line)) {
        if(line.find("MSG,3") != -1) {
            // cout << "1" << line << "\n";
            numProcess++;
            for(char ch : line) {
                cout << ch << endl;
            }
            break;
        }
    }
    MyFile.close();
}

void createSMTProg::readSegTrace(int start, int end, int eps) {
    int i, j, num = 1;
    start = max(start - eps, 0);
    // cout << start << " : " << end << endl;
    for(i = 0; i < numProcess; i++) {
        numEvent[i] = 0;
        for(j = 0; j < numEventAll[i]; j++) {
            if(ETimeAll[i][j][0] >= start && ETimeAll[i][j][0] <= end) {
                numEvent[i]++;
                EName[i][j] = to_string(num++);
                EType[i][j] = ETypeAll[i][j];
                EValue[i][j] = EValueAll[i][j];
                ETime[i][j][0] = ETimeAll[i][j][0];
                ETime[i][j][1] = ETimeAll[i][j][1];
                ETime[i][j][2] = ETimeAll[i][j][2];
            }
        }
    }
    // cout << num << "\n";
    for(int i = 0; i < numProcess; i++) {
        for(int j = 0; j < numEvent[i]; j++) {
            if(EType[i][j] == "R") {
                for(int k = 0; k < numProcess; k ++) {
                    for(int l = 0; l < numEvent[k]; l++) {
                        if(EType[k][l] == "S" && ETime[i][j][0] == stoi(EValue[i][j].substr(0, EValue[i][j].find(",")))) {
                            hbSet[numHbSet][0] = stoi(EName[i][j]);
                            hbSet[numHbSet++][1] = stoi(EName[k][l]);
                            // cout << "Happen Before: " << hbSet[numHbSet - 1][0] << " : " << hbSet[numHbSet - 1][1] << "\n";
                        }
                    }
                }
            }
            for(int k = 0; k < numProcess; k++) {
                for(int l = 0; l < numEvent[k]; l++) {
                    if(abs(ETime[i][j][0] - ETime[k][l][0]) >= eps && i != k) {
                        hbSet[numHbSet][0] = stoi(EName[i][j]);
                        hbSet[numHbSet++][1] = stoi(EName[k][l]);
                        // cout << "Happen Before: " << hbSet[numHbSet - 1][0] << " : " << hbSet[numHbSet - 1][1] << "\n";
                    }
                }
            }
        }
    }
}

void createSMTProg::makeSMTProg(int eps, int segLength, int maxTime, string phi, int numThreads)
{
    readFile();

    // string formulaList[10];
    // int numFormula = 3;
    // string currentFormula[10];

    // int resultMat[500][10];
    // for(int  i = 0; i < int (maxTime/ segLength) + 1; i++)
    // {
    //     for(int j = 0; j < 2 * numFormula; j++)
    //     {
    //             resultMat[i][j] = 0;
    //     }
    // }

    // int numSegment = -1;
    // // numThreads = 1;
    // int segmentL = maxTime / numThreads;
    // // cout << "total numThreads: " << numThreads << endl;
    // #pragma omp parallel num_threads(numThreads)
    // {
    //     #pragma omp for
    //     for(int startT = 0; startT < maxTime; startT = startT + segmentL) {
    //         // cout << "Thread: " << startT << " : " << startT + segmentL << endl;

    //         int start, end;
    //         end = startT;

    //         while(end < startT + segmentL) {
    //             numSegment++;
    //             start = end;
    //             end += segLength;
    //             // cout << start << " : " << end << "\n";

    //             readSegTrace(start, end, eps);
    //             numFormula = 1;
    //             for (int i = 0; i < 2 * numFormula; i++) {
    //                     // cout << statesList[i] << " : " << statesList[j] << "\n";
    //                     resultMat[numSegment][i] = solveSMTProg(formulaList[i % 2]);
    //             }
    //         }
    //     }
    }

    // for(int i = 0; i < numSegment; i++) {
    //     string newCurrentFormula[10];
    //     int numNewCurrentFormula = 0;
    //     for(int j = 0; j < numCurrentFormula; j++) {
    //         for(int k = 0; k < numFormula; k++) {
    //             if(formulaList[k] == currentFormula[j]) {
    //                 for(int l = 0; l < numFormula; l++) {
    //                     if(resultMat[i][k][l] == 1)
    //                         newCurrentFormula[numNewCurrentFormula ++] = formulaList[l];
    //                 }
    //             }
    //         }
    //     }
    //     numCurrentFormula = 0;
    //     for(int j = 0; j < numNewCurrentFormula; j++) {
    //         currentFormula[numCurrentFormula++] = newCurrentFormula[j];
    //     }
    // }
}

bool createSMTProg::solveSMTProg(string subFormula) {
    context c;

    solver s(c);

    expr_vector eventList(c);

    int totNumEvents = 0;
    for(int i = 0; i < numProcess; i++)
        totNumEvents += numEvent[i];

    for(int i = 0; i < pow(totNumEvents, 2); i++)
    {
        std::string str = "eventList" + std::to_string(i);
        eventList.push_back(c.bv_const(str.c_str(), totNumEvents));
        s.add(eventList[i] == i);
    }

    func_decl f = z3::function("f", c.int_sort(), c.bv_sort(totNumEvents));

    s.add(f(0) == eventList[0]);

    for (int i = 0; i <= totNumEvents; i++)
    {
        expr_vector event_range(c);

        for (int j = 0; j < pow(totNumEvents, 2); j++)
        {
            event_range.push_back(f(i) == eventList[j]);
        }

        s.add(mk_or(event_range));
    }

    for(int i = 0; i < totNumEvents; i++)
    {
        expr_vector event_order(c);

        for(int j = 1; j < pow(totNumEvents, 2); j = j * 2)
            event_order.push_back(bv2int(f(i + 1), false) - bv2int(f(i), false) == j);

        // s.add(f(i + 1) > f(i));
        s.add(mk_or(event_order));
    }

    s.add(f(totNumEvents) == eventList[(int)pow(totNumEvents, 2) - 1]);

    expr x = c.int_const("x");
    s.add(0 <= x && x <= totNumEvents);

    expr_vector b1(c);
    expr_vector b2(c);

    for(int i = 0; i < numHbSet; i++)
    {
        std::string str = "b1" + std::to_string(i);
        b1.push_back(c.bv_const(str.c_str(), totNumEvents));
        str = "b2" + std::to_string(i);
        b2.push_back(c.bv_const(str.c_str(), totNumEvents));
        int p0 = (int)pow(2, hbSet[i][0] - 1);
        int p1 = (int)pow(2, hbSet[i][1] - 1);

        s.add(b1[i] == p1);
        s.add(b2[i] == p0);
        s.add(forall(x, implies(bv2int(f(x) & b1[i], false) != 0, bv2int(f(x) & b2[i], false) != 0)));
    }

    // func_decl_vector frontier(c);

    // for(int i = 0; i < numProcess; i++)
    // {
    //     std::string str = "frontier" + std::to_string(i);
    //     frontier.push_back(z3::function(str.c_str(), c.int_sort(), c.bv_sort(totNumEvents)));

    //     expr_vector frontier_order(c);

    //     for(int k = 0; k < totNumEvents; k++)
    //     {
    //         for(int j = 0; j < numEvent[i]; j++)
    //         {
    //             frontier_order.push_back(frontier[i](k) == EName[i][j]);
    //         }
    //         s.add(mk_or(frontier_order));
    //     }
    // }

    if(s.check() == sat)
    {
//        std::string verdict = std::to_string(i) + " : Sat";
        // std::cout << "Sat" << std::endl;
        return true;

//        model m = s.get_model();
//        std::cout << m << std::endl;
    }
    else
    {
//        std::string verdict = std::to_string(i) + " : Unsat";
        // std::cout << "Unsat" << std::endl;
        return false;

    }
    Z3_reset_memory();
    s.reset();
}



class mainProg {

private:
    int numCores;
    float average(float*, int);
    float standardDeviation(float*, int, float);

public:
    void setNumCores(int);
    mainProg();
    mainProg(int);
    void execute();
};

void mainProg::setNumCores(int num) {
    numCores = num;
}

mainProg::mainProg() {
    numCores = 1;
}

mainProg::mainProg(int num) {
    numCores = num;
}

float mainProg::average(float num[], int n)
{
    float sum = 0;
    for(int i = 0; i < n; i++)
        sum += num[i];
    return sum/n;
}

float mainProg::standardDeviation(float num[], int n, float m)
{
    float sum = 0;
    for(int i = 0; i < n; i++)
    {
        sum += pow(num[i] - m, 2);
    }
    return sqrt(sum/n);
}

void mainProg::execute() {
    // cout << "hi";
//    SynthExp synthExp("trace2_20.txt", "pq", 2, 20);
//    synthExp.genTrace();
    int numIter = 1;
    float result[10];
    float meanTime, sdTime;
    clock_t runTime;

    ofstream outFile;
    outFile.open("report.csv");
    outFile << "numProcess, segLength, compLength, eps, eventRate, monTime, monTime-SD, progTime, progTime-SD" << endl;
    for(int j = 0;j < 1; j++) {
        for(int numProcess = 2; numProcess <= 2; numProcess++) {
            for(int segLength = 10; segLength <= 100; segLength = segLength + 2) {
                int compLength = 200;
                int eps = 10;
                // int numProcess = 2;
                // int segLength = 10;
                int eventRate = 10;

                string traceFile = "sbs1.txt"
                // string command = "python synth-system.py " + to_string(numProcess) + " " + to_string(compLength) + " " + traceFile + " " + to_string(eventRate);
                // system(command.c_str());

                cout << "numProcess: " << numProcess << " compLength: " << compLength << " segLength: " << segLength << " eps: " << eps << " eventRate: " << eventRate << "\n";
                outFile << numProcess << ", " << segLength << ", " << compLength << ", " << eps << ", " << eventRate << ", ";

                for(int i = 0; i < numIter; i++)
                {
                    runTime = clock();
                    createSMTMon smtMon(traceFile);
                    smtMon.makeSMTMon(eps, segLength, compLength, numCores);
                    runTime = clock() - runTime;
                    result[i] = (float) runTime / CLOCKS_PER_SEC;
                    // cout << result[i] << endl;
                }

                meanTime = average(result, numIter);
                sdTime = standardDeviation(result, numIter, meanTime);
                cout << "MonitorTime: " << meanTime << ", SD: " << sdTime << endl;
                outFile << meanTime << ", " << sdTime << ", ";

                for(int i = 0; i < numIter; i++)
                {
                    runTime = clock();
                    createSMTProg smtProg(traceFile);
                    smtProg.makeSMTProg(eps, segLength, compLength, "<> r -> ( ! p U r)", numCores);
                    runTime = clock() - runTime;
                    result[i] = (float) runTime / CLOCKS_PER_SEC;
                    // cout << result[i] << endl;
                }

                meanTime = average(result, numIter);
                sdTime = standardDeviation(result, numIter, meanTime);
                cout << "ProgressionTime: " << meanTime << ", SD: " << sdTime << endl;
                outFile << meanTime << ", " << sdTime << endl;
            }
        }
    }
    outFile.close();
}

int main() {
    cout << "Welcome!" << endl;
    mainProg obj(1);
    obj.execute();

    return 1;
}
