#include "z3++.h"
#include <string>
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <time.h>
#include <omp.h>
#include <stdlib.h>
#include <string>

using namespace std;
using namespace z3;

class createSMTProg {
private:
    string traceFileName;
    string ENameAll[150];
    float EAltAll[150][20];
    float ELatlAll[150][20];
    float ELogtAll[150][20];
    int ETimeAll[150][20];
    int numProcess, numEventAll[150];

    int EName[150][20];
    float EAlt[150][20];
    float ELatl[150][20];
    float ELogt[150][20];
    int ETime[150][20];
    int numEvent[150];

    int hbSet[500][2];
    int numHbSet;

    int maxTime;

    void readFile(string, string);
    void cleanData();

public:
    createSMTProg();
    createSMTProg(string);
    void makeSMTProg(int, int, int);
    bool solveSMTProg(string);
};

createSMTProg::createSMTProg() {
    traceFileName = "sbs1.txt";
    cleanData();
}

createSMTProg::createSMTProg(string one) {
    traceFileName = one;
    cleanData();
}

void createSMTProg::cleanData()
{
    maxTime = 0;
    int i, j;
    for(i = 0; i < 150; i++) {
        ENameAll[i] = "";
        for(j = 0; j < 20; j++) {
            EAltAll[i][j] = 0;
            ELatlAll[i][j] = 0;
            ELogtAll[i][j] = 0;
            ETimeAll[i][j] = 0;

            EName[i][j] = 0;
            EAlt[i][j] = 0;
            ELatl[i][j] = 0;
            ELogt[i][j] = 0;
            ETime[i][j] = 0;
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

void createSMTProg::readFile(string start, string end) {
    ifstream MyFile(traceFileName);
    string line;
    int proc = 0, event = -1, num = 0, wordNum;
    string word = "";
    double startTime = 0, endTime = 0, wordTime;
    for(char ch1 : start) {
        if(ch1 >= '0' && ch1 <= '9')
            startTime = startTime * 10 + (ch1 - '0');
    }
    for(char ch1 : end) {
        if(ch1 >= '0' && ch1 <= '9')
            endTime = endTime * 10 + (ch1 - '0');
    }
    endTime = endTime - startTime;

    while(getline(MyFile, line)) {
        if(line.find("MSG,3") != -1) {
            wordNum = 0;
            for(char ch : line) {
                if(ch == ','){
                    wordNum++;
                    if(wordNum == 5) {
                        // cout << "Name" << numProcess << ": " << word << endl;
                        proc = -1;
                        for(int i = 0; i < numProcess; i++) {
                            if(ENameAll[i].compare(word) == 0)
                                proc = i;
                        }
                        if(proc == -1) {
                            proc = numProcess;
                        }
                    }
                    else if(wordNum == 8) {
                        wordTime = 0;
                        for(char ch1 : word) {
                            if(ch1 >= '0' && ch1 <= '9')
                                wordTime = wordTime * 10 + (ch1 - '0');
                        }
                        wordTime = wordTime - startTime;
                        if(endTime < wordTime) {
                            break;
                        }
                        else if(wordTime < 0) {
                            continue;
                        }
                        // cout << "Time: " << word << " : " << wordTime << endl;
                        ETimeAll[proc][numEventAll[proc]++] = wordTime;
                        num++;
                        maxTime = wordTime;
                        ENameAll[proc] = word;
                        if(numEventAll[proc] == 1)
                            numProcess++;
                    }
                    else if(wordNum == 12) {
                        // cout << "Altitude: " << word << endl;
                        EAltAll[proc][numEventAll[proc] - 1] = stof(word);
                    }
                    else if(wordNum == 15) {
                        // cout << "Latitude: "<< word << endl;
                        ELatlAll[proc][numEventAll[proc] - 1] = stof(word);
                    }
                    else if(wordNum == 16) {
                        // cout << "Longitude: "<< word << endl;
                        ELogtAll[proc][numEventAll[proc] - 1] = stof(word);
                    }
                    word = "";
                }
                else
                    word = word + ch;
            }
            if(endTime < wordTime)
                break;
        }
    }
    // cout << numProcess << ":" << num << endl;
    MyFile.close();
}

void createSMTProg::makeSMTProg(int eps, int segLength, int numThreads)
{
    int h = 17, m = 26, s = 52;

    while(true) {
        string startTime = to_string(h) + ":";
        if(m == 0)
            startTime = startTime + "00:";
        else if(m < 10)
            startTime = startTime + "0" + to_string(m) + ":";
        else
            startTime = startTime + to_string(m) + ":";
        if(s < 10)
            startTime = startTime + "0" + to_string(s) + ".000";
        else
            startTime = startTime + to_string(s) + ".000";
        cout << startTime << " : ";
        s = s + 5;
        if(s >= 60){
            s = s - 60;
            m = m + 1;
            if(m >= 60){
                m = m - 60;
                h++;
            }
        }
        string endTime = to_string(h) + ":";
        if(m == 0)
            endTime = endTime + "00:";
        else if(m < 10)
            endTime = endTime + "0" + to_string(m) + ":";
        else
            endTime = endTime + to_string(m) + ":";
        if(s < 10)
            endTime = endTime + "0" + to_string(s) + ".000";
        else
            endTime = endTime + to_string(s) + ".000";
        cout << endTime << endl;
        readFile(startTime, endTime);

        int numSegment = -1;
        int segmentL = maxTime / numThreads;
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

                    int i, j, num = 1;
                    int start1 = max(start - eps, 0);
                    // cout << start << " : " << end << endl;
                    for(i = 0; i < numProcess; i++) {
                        numEvent[i] = 0;
                        for(j = 0; j < numEventAll[i]; j++) {
                            if(ETimeAll[i][j] >= start1 && ETimeAll[i][j] <= end) {
                                numEvent[i]++;
                                EName[i][numEvent[i] - 1] = num++;
                                EAlt[i][numEvent[i] - 1] = EAltAll[i][j];
                                ELogt[i][numEvent[i] - 1] = ELogtAll[i][j];
                                ELatl[i][numEvent[i] - 1] = ELatlAll[i][j];
                                ETime[i][numEvent[i] - 1] = ETimeAll[i][j];
                            }
                        }
                    }
                    // cout << num << "\n";
                    for(int i = 0; i < numProcess; i++) {
                        for(int j = 0; j < numEvent[i]; j++) {
                            for(int k = 0; k < numProcess; k++) {
                                for(int l = 0; l < numEvent[k]; l++) {
                                    if(abs(ETime[i][j] - ETime[k][l]) >= eps && i != k) {
                                        hbSet[numHbSet][0] = EName[i][j];
                                        hbSet[numHbSet++][1] = EName[k][l];
                                        // cout << "Happen Before: " << hbSet[numHbSet - 1][0] << " : " << hbSet[numHbSet - 1][1] << "\n";
                                    }
                                }
                            }
                        }
                    }

                    int numFormula = 1;
                    for (int i = 0; i < 2 * numFormula; i++) {
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
                            // cout << "Sat" << std::endl;
                            // resultMat[numSegment][i] = true;

                            //        model m = s.get_model();
                            //        std::cout << m << std::endl;
                        }
                        else
                        {
                            //        std::string verdict = std::to_string(i) + " : Unsat";
                            // cout << "Unsat" << std::endl;
                            // resultMat[numSegment][i] = false;

                        }
                        Z3_reset_memory();
                        s.reset();
                    }
                }
            }
        }

        cleanData();

        if(h == 17 && m == 26 && s == 57)
            break;
    }
}

bool createSMTProg::solveSMTProg(string subFormula) {
    return true;
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
            for(int segLength = 1000; segLength <= 1000; segLength = segLength + 2) {
                int compLength = 200;
                int eps = 1000;
                // int numProcess = 2;
                // int segLength = 10;
                int eventRate = 10;

                string traceFile = "sbs1.txt";

                // cout << "numProcess: " << numProcess << " compLength: " << compLength << " segLength: " << segLength << " eps: " << eps << " eventRate: " << eventRate << "\n";
                outFile << numProcess << ", " << segLength << ", " << compLength << ", " << eps << ", " << eventRate << ", ";

                for(int i = 0; i < numIter; i++)
                {
                    runTime = clock();
                    createSMTProg smtProg(traceFile);
                    smtProg.makeSMTProg(eps, segLength, numCores);
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
