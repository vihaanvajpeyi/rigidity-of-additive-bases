// Simulated annealing: find A subset Z_m, |A|=k, with 1 <= R_A(n) <= U for all n (ordered counts).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
static int m,k,U; static int A[4096], inA[8192], r[8192];
static unsigned long long s[2]; static double TLO=0.35, THI=1.25;
static inline unsigned long long rnd(void){ unsigned long long s1=s[0], s0=s[1]; s[0]=s0; s1^=s1<<23; s[1]=s1^s0^(s1>>17)^(s0>>26); return s[1]+s0; }
static inline double urand(void){ return (rnd()>>11)*(1.0/9007199254740992.0); }
static int W4=1; static inline int pen(int v){ return v==0 ? 3 : (v>U ? v-U : (v==4 ? W4 : 0)); }
int main(int argc,char**argv){
  m=atoi(argv[1]); k=atoi(argv[2]); U=atoi(argv[3]); double secs=atof(argv[4]); unsigned seed=argc>5?atoi(argv[5]):1; if(argc>7){TLO=atof(argv[6]);THI=atof(argv[7]);}
  s[0]=0x9E3779B97F4A7C15ULL^seed; s[1]=0xD1B54A32D192ED03ULL*(seed+1); for(int i=0;i<20;i++) rnd();
  clock_t t0=clock(); long long it=0; int best=1<<30; int restarts=0;
  while(1){
    memset(inA,0,sizeof inA); memset(r,0,sizeof r);
    int c=0; while(c<k){ int x=rnd()%m; if(!inA[x]){inA[x]=1;A[c++]=x;} }
    for(int i=0;i<k;i++) for(int j=0;j<k;j++){ int t=A[i]+A[j]; if(t>=m)t-=m; r[t]++; }
    int cost=0; for(int n=0;n<m;n++) cost+=pen(r[n]);
    double T=THI; long long lim= (long long)m*k*400; 
    for(long long q=0;q<lim && cost>0;q++){
      it++;
      int i=rnd()%k; int a=A[i]; int b=rnd()%m; if(inA[b]) continue;
      // remove a: sums a+x for x in A (x!=a) lose 2, 2a loses 1 ; add b: b+x (x in A\{a}) gain 2, 2b gains 1
      int d=0;
      for(int j=0;j<k;j++){ if(j==i) continue; int t=a+A[j]; if(t>=m)t-=m; d-=pen(r[t]); r[t]-=2; d+=pen(r[t]); }
      { int t=(2*a)%m; d-=pen(r[t]); r[t]-=1; d+=pen(r[t]); }
      for(int j=0;j<k;j++){ if(j==i) continue; int t=b+A[j]; if(t>=m)t-=m; d-=pen(r[t]); r[t]+=2; d+=pen(r[t]); }
      { int t=(2*b)%m; d-=pen(r[t]); r[t]+=1; d+=pen(r[t]); }
      if(d<=0 || urand()<exp(-d/T)){ inA[a]=0; inA[b]=1; A[i]=b; cost+=d; if(cost<best) best=cost; }
      else { // undo
        for(int j=0;j<k;j++){ if(j==i) continue; int t=b+A[j]; if(t>=m)t-=m; r[t]-=2; }
        r[(2*b)%m]-=1;
        for(int j=0;j<k;j++){ if(j==i) continue; int t=a+A[j]; if(t>=m)t-=m; r[t]+=2; }
        r[(2*a)%m]+=1;
      }
      if((q&0xFFFF)==0){ T = TLO + (THI-TLO)*urand(); if((double)(clock()-t0)/CLOCKS_PER_SEC>secs) goto done; }
    }
    if(cost==0){ 
      // verify from scratch
      int rr[8192]; memset(rr,0,sizeof rr); for(int i=0;i<k;i++) for(int j=0;j<k;j++) rr[(A[i]+A[j])%m]++;
      int mn=1<<30,mx=0; for(int n=0;n<m;n++){ if(rr[n]<mn)mn=rr[n]; if(rr[n]>mx)mx=rr[n]; }
      // sort
      for(int i=0;i<k;i++) for(int j=i+1;j<k;j++) if(A[j]<A[i]){int t=A[i];A[i]=A[j];A[j]=t;}
      printf("FOUND m=%d k=%d U=%d min=%d max=%d time=%.1f A=[",m,k,U,mn,mx,(double)(clock()-t0)/CLOCKS_PER_SEC);
      for(int i=0;i<k;i++) printf("%d%s",A[i],i+1<k?",":""); printf("]\n"); return 0; }
    restarts++;
    if((double)(clock()-t0)/CLOCKS_PER_SEC>secs) break;
  }
done:
  printf("NOTFOUND m=%d k=%d U=%d best_cost=%d iters=%lld restarts=%d\n",m,k,U,best,it,restarts); return 1;
}
