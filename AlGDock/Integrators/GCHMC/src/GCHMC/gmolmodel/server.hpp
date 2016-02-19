#ifndef SERVER_H_
#define SERVER_H_


#ifndef CLIENT_FLAG
#define CLIENT_FLAG  1.0
#endif

#ifndef SERVER_FLAG
#define SERVER_FLAG  0.0
#endif

#ifndef CNT_START
#define CNT_START  1945
#endif

#ifndef CNT_STOP
#define CNT_STOP  1947
#endif

//TARGET_TYPE **indexMap;
//TARGET_TYPE *PrmToAx_po;
//TARGET_TYPE *MMTkToPrm_po;
#ifndef AxToPrm
#define AxToPrm(i) indexMap[i][1]
#endif

#ifndef PrmToAx
#define PrmToAx(i) PrmToAx_po[i]
#endif
#ifndef PrmToMMTk
#define PrmToMMTk(i) indexMap[i][2]
#endif

#ifndef MMTkToPrm
#define MMTkToPrm(i) MMTkToPrm_po[i]
#endif
//#define AxToMMTk(i) PrmToMMTk(AxToPrm(i))
//#define MMTkToAx(i) PrmToAx(MMTkToPrm(i))
#ifndef AxToMMTk
#define AxToMMTk(i) indexMap[(int)(indexMap[i][1])][2]
#endif
#ifndef MMTkToAx
#define MMTkToAx(i) PrmToAx_po[(int)(MMTkToPrm_po[i])]
#endif


#endif //SERVER_H_
