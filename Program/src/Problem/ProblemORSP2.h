// *******************************************************************
//      file with specific functions to solve a ORSP
// *******************************************************************
#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "../Main/Data.h"

#include <math.h>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <list>
#include <algorithm>    


// Variables declared in main.cpp
extern int n;                               // size of cromossoms

//------ DEFINITION OF TYPES OF PROBLEM SPECIFIC --------

#define ROOMTYPE_BED 0
#define ROOMTYPE_OR 1
#define ROOMTYPE_ICU 2
#define ROOMTYPE_WARD 3
#define NUM_ROOMTYPES 4

#define PERSONTYPE_STAFF 0
#define PERSONTYPE_SURGEON 1
#define PERSONTYPE_PACIENT 2
#define NUM_PERSONTYPES 3

#define EQUIPTYPE_MACHINEA 0
#define EQUIPTYPE_MACHINEB 1
#define EQUIPTYPE_MACHINEC 2
#define EQUIPTYPE_MACHINED 3
#define EQUIPTYPE_BLADES 4
#define NUM_EQUIPTYPES 5


#define DAYS_PER_WEEK 7
#define MIN_PER_DAY 1440
#define MIN_PER_WEEK 10080
#define BUSINESS_MIN_PER_DAY 540


//Static Structures

struct Equipament				// struct with event info
{
	const int id;
	const int type;             // machine, blades
    int preTime;
    int postTime;
    int numAvailable;
    Equipament(int _id, int _tp) : id(_id), type(_tp) {}
};

struct Room				// struct with room info
{
	int id;
	int type;           // bed, OR, ICU, ward
    int timeAvaiable;
    Room(int _id, int _tp) : 
        id(_id), type(_tp), timeAvaiable(0) {}
    Room(int _id, int _tp, int _ta) : 
        id(_id), type(_tp), timeAvaiable(_ta) {}
};

struct Person         // struct with person info
{
	int id;
	int type;         //0 Staff 1 Surgeon; 2 Client
    int timeAvaiable;
    Person(int _id, int _tp) :
        id(_id), type(_tp), timeAvaiable(0) {}
    Person(int _id, int _tp, int _ta) :
        id(_id), type(_tp), timeAvaiable(_ta) {}
};

struct Event			// struct with event info
{
	const int id;
    const int jobID;
    const int uniqID;
    const int roomType;
    int prepTime;
    int duration;
    int blockingLimit;
    int transferTime;
    int postTime;
    bool businessHours;
    int equipaments[NUM_EQUIPTYPES];
    int personal[NUM_PERSONTYPES];

    Event(int _id, int _jID, int _uID, int _rt) : 
     id(_id), jobID(_jID), uniqID(_uID), roomType(_rt), prepTime(0),
     duration(0), transferTime(0), postTime(0), businessHours(false) {
        for (int i=NUM_EQUIPTYPES-1;i>=0;--i) equipaments[i]=0;
        for (int i=NUM_PERSONTYPES-1;i>=0;--i) personal[i]=0;
    }
};

struct Surgery		// struct with node informations
{
	const int id;
	const int clientID;
    const int surgeonID;
    std::vector <Event> events;

    Surgery(int _id, int _clt, int _sur) : 
         id(_id),  clientID(_clt), surgeonID(_sur) {}
};

//Dynamic Structures

struct RoomReserve	// struct with room reservation informations
{
    Event *event;
    int reserveTime;
    int startTime;
    int endTime;
    int blocking;
    std::list <Person*> people;
};

struct RoomSchedule		// struct with room schedule
{
    Room *room;
    std::list <RoomReserve*> schedule;
    int timeAvaiable;
};

struct PersonSchedule  // struct with person schedule
{
    Person *person;
    std::list <RoomReserve*> schedule;
    int timeAvaiable;
};

struct EquipamentSchedule // struct for multiple equipaments allocation
{
    Equipament *equip;
    std::vector <RoomReserve*> alloc;
};

struct InsertionData // struct for 
{
    std::vector<int> timesReserve;
    std::vector<int> timesStart;
    std::vector<int> timesEnd;
    std::vector<int> blocking;
    std::vector<RoomSchedule*> lstRs;
    std::vector< std::vector<PersonSchedule*> > lstPs;
    InsertionData(Surgery &surg) {
        const int nEvents = surg.events.size();
        timesReserve.resize(nEvents);
        timesStart.resize(nEvents);
        timesEnd.resize(nEvents);
        blocking.resize(nEvents);
        lstRs.resize(nEvents);
        lstPs.resize(nEvents);
        for (int i = 0; i < nEvents; i++) {
            int nPer = surg.events[i].personal[ PERSONTYPE_STAFF ] + 1;
            if (surg.events[i].roomType == ROOMTYPE_OR) nPer++;
            lstPs[i].resize(nPer);
        }
    }
};

struct DecodedSolution
{
    std::vector <RoomSchedule*> beds;
    std::vector <RoomSchedule*> opRms;
    std::vector <RoomSchedule*> icus;
    std::vector <RoomSchedule*> wards;
    std::vector <PersonSchedule*> nurses;
    std::vector <PersonSchedule*> surgeons;
    std::vector <PersonSchedule*> pacients;
    std::vector <EquipamentSchedule> equipaments;
    std::vector <RoomSchedule> AroomS;
    std::vector <PersonSchedule> ApersonS;
    std::vector <RoomReserve> ARreservs;
    int numSurgeries;
    int totalTime;
};


int planingHorizon;
int totEvents;
bool generateString;
int minStartTime;

//-------------------------- FUNCTIONS OF SPECIFIC PROBLEM --------------------------

//------ DEFINITION OF GLOBAL CONSTANTS AND VARIABLES OF SPECIFIC PROBLEM  --------

static std::vector <Equipament> equipaments;
static std::vector <Person> nurses;
static std::vector <Person> surgeons;
static std::vector <Person> pacients;
static std::vector <Room> beds;
static std::vector <Room> opRms;
static std::vector <Room> icus;
static std::vector <Room> wards;
static std::vector <Surgery> surgeries;
static std::vector <Event> events;

//-------------------------- FUNCTIONS OF SPECIFIC PROBLEM --------------------------

/************************************************************************************
 Method: insert_in_order
 Description: Inserts a element on a sorted list
*************************************************************************************/
void insert_in_order(std::list<RoomReserve*> &my_list, RoomReserve *element);

/************************************************************************************
 Method: CalcReqNumber
 Description: Calculate the minimum amount of equipaments to be allocated
*************************************************************************************/
int CalcReqNumber(EquipamentSchedule &equipS);

/************************************************************************************
 Method: GreedInsertion
 Description: Greedly inserts a client at the earliest avaiable time slots
*************************************************************************************/
bool GreedInsertion(DecodedSolution &dSol, Surgery &surg, const int srType);

/************************************************************************************
 Method: FreeMemoryProblem
 Description: Free local memory allocate by Problem
*************************************************************************************/
void FreeMemoryProblem();


std::vector< std::string > splitString(const char* input, const char sep) {
    std::stringstream test(input);
    std::string segment;
    std::vector<std::string> seglist;
    while(std::getline(test, segment, sep)) seglist.push_back(segment);
    return seglist;
}

void printPersonSchedule(const PersonSchedule *ps, std::string &strOut) {
    char buffer[500];
    sprintf(buffer, "Person %d-%d\t", ps->person->type+1, ps->person->id+1); strOut.append(buffer);
    for (RoomReserve *rr : ps->schedule) {
        sprintf(buffer, "%d-%d %d-%d\t", rr->event->jobID+1, rr->event->id+1, rr->startTime, rr->endTime);
        strOut.append(buffer);
    }
    strOut.append("\n");
}

void printRoomSchedule(const RoomSchedule *rs, std::string &strOut) {
    char buffer[500];
    sprintf(buffer, "Room %d-%d\t", rs->room->type+1, rs->room->id+1); strOut.append(buffer);
    for (RoomReserve *rr : rs->schedule) {
        sprintf(buffer, "%d-%d %d-%d-%d\t", rr->event->jobID+1, rr->event->id+1, rr->startTime, rr->endTime, rr->blocking);
        strOut.append(buffer);
    }
    strOut.append("\n");
}

void printSolutionData(const DecodedSolution &dSol, std::string &strOut) {
    if (generateString == false)  return;
    printf("Here\n");
    char buffer[500]; strOut.clear();
    sprintf(buffer, "makespan:%d\tnumSurgeries:%d\n", dSol.totalTime, dSol.numSurgeries); strOut.append(buffer);
    for (const RoomSchedule *rs : dSol.beds) printRoomSchedule(rs, strOut);
    for (const RoomSchedule *rs : dSol.opRms) printRoomSchedule(rs, strOut);
    for (const RoomSchedule *rs : dSol.icus) printRoomSchedule(rs, strOut);
    for (const RoomSchedule *rs : dSol.wards) printRoomSchedule(rs, strOut);
}

void printInputData() {
    for (Surgery &sur : surgeries) {
        printf("jbID:%d clientID:%d surgeonID:%d nEvents:%d\n", sur.id, sur.clientID, sur.surgeonID, (int)sur.events.size());
        for (Event &evt : sur.events) {
            printf("eventID:%d roomType:%d prTime:%d duration:%d trasfTime:%d pstTime:%d BH:%d\n", evt.id, evt.roomType, 
                   evt.prepTime,evt.duration,evt.transferTime, evt.postTime, (int)evt.businessHours);
        }
    }
    printf("\n\n");
    getchar();
}

Surgery* getAddSurgery(const int _Surgeryid, const int _pacientId, const int _SurgeonId) {
    for (Surgery &job : surgeries) {
        if (job.id == _Surgeryid) return &job; 
    }
    /*Person *psur = NULL;
    for (Person &pSurg : surgeons) { if (pSurg.id == _SurgeonId) { psur = &pSurg; break; } }
    if (psur == NULL) {
        surgeons.push_back(Person(_SurgeonId, PERSONTYPE_SURGEON));
        psur = &surgeons.back();
    }
    Person *pPac = NULL;
    for (Person &bClient : pacients) { if (bClient.id == _pacientId) { pPac = &bClient; break; } }
    if (pPac ==NULL) {
        pacients.push_back(Person(_pacientId, PERSONTYPE_PACIENT));
        pPac = &pacients.back();
    }*/
    if ((int)surgeons.size() <= _Surgeryid) surgeons.push_back(Person(_SurgeonId, PERSONTYPE_SURGEON));
    if ((int)pacients.size() <= _Surgeryid) pacients.push_back(Person(_pacientId, PERSONTYPE_PACIENT));
    //printAllData();
    surgeries.push_back(Surgery( _Surgeryid, _pacientId, _SurgeonId));
    Surgery *sur = &surgeries.back();
    //printf("%d %d-%d %d-%d\n", sur->id, pPac->id, sur->client.id, psur->id, sur->surgeon.id);
    return sur;
}

Event* getAddEvent(const int _Surgeryid, const int _EventId, const int roomType) {
    Surgery *sur = NULL;
    for (Surgery &job : surgeries) {
        if (job.id == _Surgeryid) {sur = &job; break;}
    }
    for (Event &bEve : sur->events) {
        if (bEve.id == _EventId) { return &bEve; }
    }
    sur->events.push_back( Event(_EventId, _Surgeryid, totEvents++, roomType) );
    return &sur->events.back();
}

void addRoom(const int _RType, const int _RoomId) {
    std::vector <Room> *vecRm = NULL;
    if (_RType == ROOMTYPE_BED ) vecRm = &beds;
    if (_RType == ROOMTYPE_OR ) vecRm = &opRms;
    if (_RType == ROOMTYPE_ICU ) vecRm = &icus;
    if (_RType == ROOMTYPE_WARD ) vecRm = &wards;
    for (Room &rm : *vecRm) if (rm.id == _RoomId) return;
    vecRm->push_back( Room(_RoomId, _RType) );
}

void ReadData(char nameTable[]) //Australian
{
    char name[200] = "../Instances/";
    strcat(name, nameTable);

    FILE *arq;
    arq = fopen(name, "r");

    if (arq == NULL)
    {
        printf("\nERROR: File (%s) not found!\n", name);
        getchar();
        exit(1);
    }

    // => read data
    char temp[500], *r = NULL;
    totEvents = 0;
    // numDecodes = 0;
    generateString = false;
    while (!feof(arq))
    {
        r = fgets(temp, sizeof(temp), arq);
        if (r == NULL) { printf("\nERROR READING\n");  }
        //printf("%s", temp);
        std::vector< std::string > lineList = splitString(temp, ',');
        if (temp[0] == 'J') {
            std::vector< std::string > jbAndEvent = splitString(lineList[1].substr(1).c_str(), '_');
            const int jobID = std::stoi( jbAndEvent[0] ) - 1;
            const int eventID = std::stoi( jbAndEvent[1] ) - 1;
            const int blkMax = std::stoi( lineList[3] );
            const int roomType = std::stoi( lineList[4].substr(1,1) ) - 1;
            const int eventDuration = std::stoi( lineList[5] );
            const int numRooms = (lineList.size()-4)/2;
            for (int i = 0; i < numRooms; i++) addRoom(roomType, i);
            getAddSurgery(jobID, jobID, jobID);
            //printf("jbID:%d evntID:%d rmTypID:%d time:%d numRooms:%d\n", jobID, eventID, roomType, eventDuration, numRooms);
            
            Event *event = getAddEvent(jobID, eventID, roomType);
            event->duration = eventDuration;
            event->blockingLimit = blkMax;
            
        }
        if (temp[0] == 'O') {
            std::vector< std::string > jbAndEvent = splitString(lineList[1].substr(1).c_str(), '_');
            const int jobID = std::stoi( jbAndEvent[0] ) - 1;
            const int eventID = std::stoi( jbAndEvent[1] ) - 1;
            const int transfDuration = std::stoi( lineList[2] );
            Event *event = getAddEvent(jobID, eventID, 0);
            event->transferTime = transfDuration;
        }
        if (temp[0] == 'P') {
            std::vector< std::string > jbAndEvent = splitString(lineList[1].substr(1).c_str(), '_');
            const int jobID = std::stoi( jbAndEvent[0] ) - 1;
            const int eventID = std::stoi( jbAndEvent[1] ) - 1;
            const int postDuration = std::stoi( lineList[2] );
            Event *event = getAddEvent(jobID, eventID, 0);
            event->postTime = postDuration;
        }
        if (temp[0] == 'B') {
            std::vector< std::string > jbAndEvent = splitString(lineList[1].substr(1).c_str(), '_');
            const int jobID = std::stoi( jbAndEvent[0] ) - 1;
            const int eventID = std::stoi( jbAndEvent[1] ) - 1;
            Event *event = getAddEvent(jobID, eventID, 0);
            event->businessHours = true;
        }
    }
    n = surgeries.size();
    planingHorizon = 7*24*60;
    //printf("\n%lu %lu %lu\n", surgeons.size(), pacients.size(), surgeries.size());
    //printf("%lu %lu %lu %lu\n", beds.size(), operationRooms.size(), icus.size(), wards.size());
    //printInputData();
    //exit(0);
}

int CalculateFitness(TSol &s) // calculate objective function
{
    int result = 0;
    /*for (RoomSchedule &rmS : s.decSol.wards) {
        if (rmS.timeAvaiable > result) result = rmS.timeAvaiable;
    }
    s.decSol.totalTime = result;
    s.ofv = result;*/
    return result;
}

// void insert_in_order(std::list<RoomReserve> &my_list, const RoomReserve &element)
void insert_in_order(std::list<RoomReserve*> &my_list, RoomReserve *element)
{
    std::list<RoomReserve*>::iterator begin = my_list.begin();
    const std::list<RoomReserve*>::iterator end = my_list.end();
    while ((begin != end) && (*begin)->endTime < element->endTime) ++begin;
    my_list.insert(begin, element);
}

// Sort TSol by random-keys
bool sortByStartTime(const RoomReserve *lhs, const RoomReserve *rhs) { return lhs->startTime < rhs->startTime; }

// int CalcReqNumber(EquipamentSchedule equipS);
int CalcReqNumber(EquipamentSchedule &equipS)
{
    std::list<RoomReserve*> buffList;
    std::vector<RoomReserve*> criticalElements;
    int maxUsage = 0, inUse = 0;
    const int equipType = equipS.equip->type;
    sort(equipS.alloc.begin(), equipS.alloc.end()-1, sortByStartTime);
    for (RoomReserve *res : equipS.alloc)
    {
        if (buffList.size() > 0)
        {
            // Remove Finished Jobs
            const int startT = res->endTime;
            do {
                RoomReserve *frt = buffList.front();
                if (frt->endTime < startT) {
                    inUse -= frt->event->equipaments[equipType];
                    buffList.pop_front();
                }
                else break;
            } while (true);
        }
        insert_in_order(buffList, res);
        inUse += res->event->equipaments[equipType];
        if (maxUsage < inUse)
        {
            // Store the moment of max utilization
            maxUsage = inUse;
            std::copy(buffList.begin(), buffList.end(), criticalElements.begin());
        }
    }
    return maxUsage;
}

void startPersonSchedule(PersonSchedule *ps, Person &per) {
    ps->person = & per; ps->timeAvaiable = per.timeAvaiable;  ps->schedule.clear();
}

void startRoomSchedule(RoomSchedule *rs, Room &rm) {
    rs->room = & rm; rs->timeAvaiable = rm.timeAvaiable; rs->schedule.clear();
}

void getEmptySolution(DecodedSolution &newSol) {
    newSol.ARreservs.resize(totEvents); int crrID = totEvents;
    newSol.AroomS.resize(beds.size() + opRms.size() + icus.size() + wards.size()); int contArs = 0;
    newSol.ApersonS.resize(nurses.size() + pacients.size() + surgeons.size()); int contAps = 0;

    newSol.numSurgeries = 0; newSol.totalTime = 0;
    newSol.nurses.resize(nurses.size()); newSol.pacients.resize(pacients.size()); newSol.surgeons.resize(surgeons.size());
    newSol.beds.resize(beds.size()); newSol.opRms.resize(opRms.size()); 
    newSol.icus.resize(icus.size()); newSol.wards.resize(wards.size());
    for (int i = nurses.size()-1; i >= 0; --i) {
        newSol.nurses[i] = &newSol.ApersonS[contAps++]; 
        startPersonSchedule(newSol.nurses[i], nurses[i]);
    }
    for (int i = pacients.size()-1; i >= 0; --i) {
        newSol.pacients[i] = &newSol.ApersonS[contAps++]; 
        startPersonSchedule(newSol.pacients[i], pacients[i]);
    }
    for (int i = surgeons.size()-1; i >= 0; --i) {
        newSol.surgeons[i] = &newSol.ApersonS[contAps++]; 
        startPersonSchedule(newSol.surgeons[i], surgeons[i]);
    }
    for (int i = beds.size()-1; i >= 0; --i) {
        newSol.beds[i] = &newSol.AroomS[contArs++]; 
        startRoomSchedule(newSol.beds[i], beds[i]);
    }
    for (int i = opRms.size()-1; i >= 0; --i) {
        newSol.opRms[i] = &newSol.AroomS[contArs++]; 
        startRoomSchedule(newSol.opRms[i], opRms[i]);
    }
    for (int i = icus.size()-1; i >= 0; --i) {
        newSol.icus[i] = &newSol.AroomS[contArs++]; 
        startRoomSchedule(newSol.icus[i], icus[i]);
    }    
    for (int i = wards.size()-1; i >= 0; --i) {
        newSol.wards[i] = &newSol.AroomS[contArs++]; 
        startRoomSchedule(newSol.wards[i], wards[i]);
    }
    for (int i = surgeries.size()-1; i >= 0; --i)  for (int j = surgeries[i].events.size()-1; j>=0; --j) {
        newSol.ARreservs[--crrID].event = & surgeries[i].events[j];
    }
}

void addSurgery(InsertionData &insData, DecodedSolution &decSol, Surgery &surg)
{
    const int nEvents = surg.events.size();
    for (int i = 0; i < nEvents; i++) {
        Event &evt = surg.events[i];
        RoomReserve *rRes = &decSol.ARreservs[evt.uniqID];
        rRes->startTime = insData.timesStart[i];
        rRes->endTime   = insData.timesEnd[i];
        rRes->blocking  = insData.blocking[i];
        insData.lstRs[i]->schedule.push_back(rRes);
        insData.lstRs[i]->timeAvaiable = fmax(insData.lstRs[i]->timeAvaiable, insData.timesEnd[i] + evt.postTime);
        for (PersonSchedule *ps : insData.lstPs[i]) {
            ps->schedule.push_back(rRes);
            ps->timeAvaiable = fmax(ps->timeAvaiable, insData.timesEnd[i] + evt.transferTime);
        }
    }
    decSol.numSurgeries++;
    decSol.totalTime = fmax(decSol.totalTime, insData.timesEnd.back());
}

// Sort the Room Schedules by a target avaiable time
bool sortRSByDelay(RoomSchedule *lhs, RoomSchedule *rhs) {
    return (fabs(lhs->timeAvaiable - minStartTime) < fabs(rhs->timeAvaiable - minStartTime)); 
}

// Sort the Person Schedules by a target avaiable time
bool sortPSByDelay(PersonSchedule *lhs, PersonSchedule *rhs) {
    return (fabs(lhs->timeAvaiable - minStartTime) < fabs(rhs->timeAvaiable - minStartTime)); 
}

RoomSchedule* pickDetRoomByType(std::vector <RoomSchedule*> *vecRS, const int typePk) {
    RoomSchedule *resp = NULL; int pickVal = 999999;
    for (RoomSchedule *rmS : *vecRS) { //Picks the room that minimizes 0-delay 1-timeAvaiable
        const int nDiff = (typePk==0) ? fabs(rmS->timeAvaiable - minStartTime) : rmS->timeAvaiable;
        //const int nDiff = rmS->timeAvaiable;
        if ( nDiff < pickVal ) { pickVal = nDiff; resp = rmS; }
    }
    return resp;
}

RoomSchedule* pickNonDetRoomByType(std::vector <RoomSchedule*> *vecRS, const int typePk) {
    int maxVal = 0; const int nElem = vecRS->size();
    // std::vector<TVecSol> vec(nElem);

    std::vector<int> sol(nElem);
    std::vector<double> rk(nElem);

    for (int i = nElem-1; i>=0; --i) { 
        sol[i] = i; 
        rk[i] = (typePk==0) ? fabs((*vecRS)[i]->timeAvaiable - minStartTime) : (*vecRS)[i]->timeAvaiable;
        maxVal = fmax(maxVal, rk[i]);
    }
    if (maxVal == 0) return (*vecRS)[ rand() % nElem ];
    //for (int i = 0; i < nElem; i++) { printf("%d-%.1lf ", vec[i].sol, vec[i].rk); } printf("\n");
    
    // sort(vec.begin(), vec.end(), sortByRk); // sort rooms by delay
    std::sort(sol.begin(), sol.end(), [&rk](int i1, int i2) {
        return rk[i1] < rk[i2];
    });


    //for (int i = 0; i < nElem; i++) { printf("%d-%.1lf ", vec[i].sol, vec[i].rk); } printf("\n");
    int totalSum = 0, maxK = fmin( nElem, 3 );
    for (int i = 0; i<maxK; i++) { 
        rk[i] = maxVal - rk[i]; 
        totalSum += rk[i];
    }
    //for (int i = 0; i < nElem; i++) { printf("%d-%.1lf ", vec[i].sol, vec[i].rk); } printf("tot: %d\n", totalSum);
    int i = 0, rng = rand() % (totalSum + 1);
    for (int acc = rk[0]; acc < rng; acc += rk[i++]);
    //printf("RNG: %d Picked: %d-%d\n", rng, i, vec[i].sol);
    return (*vecRS)[ sol[i] ];
}

bool GreedInsertion(DecodedSolution &dSol, Surgery &surg, const int srType)
{
    InsertionData insData(surg);
    minStartTime = pacients[surg.clientID].timeAvaiable;
    const int nEvents = surg.events.size();
    RoomSchedule *rsw = NULL; //int numLoops = 0;
    // printf("nEvnts %d ClientID %d surgeonID %d\n", (int)insData.lstPs.size(), surg.clientID, surg.surgeonID);
    for (bool feasInsert = false; !feasInsert;)
    {
        feasInsert = true;
        for (int i = 0; i < nEvents; i++)
        {
            Event &evt = surg.events[i];
            int remStaff = insData.lstPs[i].size();
            std::vector<RoomSchedule *> *vecRS = NULL;
            if (evt.roomType == ROOMTYPE_BED) vecRS = &dSol.beds;
            if (evt.roomType == ROOMTYPE_OR)  vecRS = &dSol.opRms;
            if (evt.roomType == ROOMTYPE_ICU) vecRS = &dSol.icus;
            if (evt.roomType == ROOMTYPE_WARD) vecRS = &dSol.wards;
            insData.lstPs[i][--remStaff] = dSol.pacients[surg.clientID];
            if (evt.roomType == ROOMTYPE_OR) insData.lstPs[i][--remStaff] = dSol.surgeons[surg.surgeonID];
            insData.lstRs[i] = (evt.roomType != ROOMTYPE_WARD) ? pickDetRoomByType(vecRS, srType) : rsw;
            //numLoops++;
            // printf("ReqStaff %d\n", remStaff);
            if (remStaff > 0)
            {
                sort(dSol.nurses.begin(), dSol.nurses.end(), sortPSByDelay); // sort staff by delay
                for (PersonSchedule *perS : dSol.nurses)
                {
                    insData.lstPs[i][--remStaff] = perS;
                    if (remStaff == 0) break;
                }
            }
            int pickVal = fmax(minStartTime, insData.lstRs[i]->timeAvaiable);
            for (PersonSchedule *perS : insData.lstPs[i])
            {
                if (perS->person->type == PERSONTYPE_STAFF && perS->timeAvaiable > pickVal) pickVal = perS->timeAvaiable;
            }
            if (evt.roomType == ROOMTYPE_BED) {
                rsw = pickDetRoomByType(&dSol.wards, srType);
                pickVal = fmax(pickVal, rsw->timeAvaiable);
            }
            insData.timesStart[i] = (evt.roomType != ROOMTYPE_WARD) ? pickVal : insData.timesStart[0];
            insData.timesStart[i] = pickVal;
            insData.timesEnd[i] = pickVal + evt.prepTime + evt.duration;
            insData.blocking[i] = 0;
            if (i > 0) {
                insData.blocking[i - 1] = pickVal - minStartTime;
                if (insData.blocking[i - 1] > surg.events[i - 1].blockingLimit) {
                    feasInsert = false;
                    minStartTime = insData.timesStart[0] + insData.blocking[i - 1];
                    break;
                }
            }
            if (evt.businessHours) {
                const int oldStartTime = insData.timesStart[i];
                const int dayWeek = (oldStartTime / MIN_PER_DAY) % DAYS_PER_WEEK;
                if (dayWeek >= 5) {
                    feasInsert = false;
                    minStartTime = ceil( (float)oldStartTime / MIN_PER_WEEK ) * MIN_PER_WEEK;
                    break;
                }
                const int startMinDay = oldStartTime % MIN_PER_DAY;
                const int endMinDay = insData.timesEnd[i] % MIN_PER_DAY;
                if (startMinDay >= 540  || endMinDay >= 540) {
                    feasInsert = false;
                    minStartTime = ceil( (float)oldStartTime / MIN_PER_DAY ) * MIN_PER_DAY;
                    break;
                }
            }
            minStartTime = insData.timesEnd[i] + evt.transferTime;
            // printf("Event %d-%d-%d RmTp:%d RoomID:%d tStart:%d tEnd:%d stfSize:%d\n", i, surg.events[i].id, surg.events[i].uniqID, surg.events[i].roomType,
            //     brmS->room->id, insData.timesStart[i], insData.timesEnd[i], (int)insData.lstPs[i].size());
        }
    }
    addSurgery(insData, dSol, surg);
    // printf("TotTime: %d numSurg: %d\n", dSol.totalTime, dSol.numSurgeries);
    //printf("numLoops: %d\n", numLoops);
    return true;
}

double Decoder(TSol s)
{
    // create a initial solution of the problem
    DecodedSolution dSol;
    getEmptySolution(dSol);

    std::vector <int> sol(n);  
    for (int j = 0; j < n; j++){ sol[j] = j;}

    // for (int i = 0; i<n; i++) {printf("%.3lf-%02d ", s.vec[i].rk, s.vec[i].sol);} printf("\n"); 
    
    s.ofv = 0;
    
    // sort the problem vector based on the values in the rk vector
    std::sort(sol.begin(), sol.end(), [&s](int i1, int i2) {
        return s.rk[i1] < s.rk[i2];
    });

    // sort(s.vec.begin(), s.vec.end(), sortByRk); // sort random-key vector 

    // for (int i = 0; i<n; i++) {printf("%.3lf-%02d ", s.vec[i].rk, s.vec[i].sol);} printf("\n"); 
    // printf("\nSize: %d", s.vec.size());
    // getchar();

    for (int j = 0; j < n && GreedInsertion(dSol, surgeries[ sol[j] ], 1); j++);      //Max Efic. per room = 0; //Pick the first avaiable room = 1
    s.ofv = (double)dSol.totalTime / MIN_PER_DAY;

    // std::string strOut;
    // printSolutionData(dSol, strOut); //s.strOut

    // print the solution in the screen
    // if (debug && print)
    // {
    //     for (int i=0; i<n; i++)
	// 	    printf("%d ", sol[i]);

    //     char buffer[500]; 
    //     std::string strOut;
    //     strOut.clear();
    //     sprintf(buffer, "\nmakespan:%d\tnumSurgeries:%d\n", dSol.totalTime, dSol.numSurgeries); strOut.append(buffer);
    //     for (const RoomSchedule *rs : dSol.beds) printRoomSchedule(rs, strOut);
    //     for (const RoomSchedule *rs : dSol.opRms) printRoomSchedule(rs, strOut);
    //     for (const RoomSchedule *rs : dSol.icus) printRoomSchedule(rs, strOut);
    //     for (const RoomSchedule *rs : dSol.wards) printRoomSchedule(rs, strOut);

    //     // Imprimir strOut na tela
    //     std::cout << strOut << std::endl;
    // }

    // // print the solution in a file
    // if (!debug && print)
    // {
    //     for (int i=0; i<n; i++)
	// 	    fprintf(arqSol,"%d ", sol[i]);

    //     char buffer[500]; 
    //     std::string strOut;
    //     strOut.clear();
    //     sprintf(buffer, "\nmakespan:%d\tnumSurgeries:%d\n", dSol.totalTime, dSol.numSurgeries); strOut.append(buffer);
    //     for (const RoomSchedule *rs : dSol.beds) printRoomSchedule(rs, strOut);
    //     for (const RoomSchedule *rs : dSol.opRms) printRoomSchedule(rs, strOut);
    //     for (const RoomSchedule *rs : dSol.icus) printRoomSchedule(rs, strOut);
    //     for (const RoomSchedule *rs : dSol.wards) printRoomSchedule(rs, strOut);
    //     fprintf(arqSol, "%s", strOut.c_str());
    // }

    return s.ofv;
}


void FreeMemoryProblem()
{
    // specific problem
    surgeries.clear(); equipaments.clear(); events.clear();
    nurses.clear(); surgeons.clear(); pacients.clear();
    beds.clear(); opRms.clear(); icus.clear(); wards.clear();
}

#endif