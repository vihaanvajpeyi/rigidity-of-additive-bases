#include <stdio.h>
#include <stdlib.h>
int m,B,R[256],sel[256],k; long long cnt; int sizes[256];
int stopfirst;
static int addel(int a){ // add a, update R; return 0 if overflow
  int ok=1;
  for(int i=0;i<k;i++){int s=(a+sel[i])%m; R[s]+=2; if(R[s]>B) ok=0;}
  int s=(2*a)%m; R[s]+=1; if(R[s]>B) ok=0;
  sel[k++]=a; return ok;}
static void remel(void){ int a=sel[--k];
  for(int i=0;i<k;i++){int s=(a+sel[i])%m; R[s]-=2;}
  R[(2*a)%m]-=1;}
int found=0;
static void rec(int next){
  // check complete
  int all=1; for(int g=0;g<m;g++) if(!R[g]){all=0;break;}
  if(all){cnt++; sizes[k]++; if(stopfirst){found=1;} }
  if(found) return;
  // pruning: k^2 <= B*m
  if((k+1)*(k+1)>B*m) return;
  for(int a=next;a<m;a++){ if(addel(a)) rec(a+1); remel(); if(found) return;}
}
int main(int c,char**v){m=atoi(v[1]);B=atoi(v[2]);stopfirst=c>3;
  for(int g=0;g<m;g++)R[g]=0; k=0; addel(0); rec(1);
  printf("m=%d B=%d count=%lld",m,B,cnt); for(int s=0;s<64;s++) if(sizes[s]) printf(" size%d:%d",s,sizes[s]); printf("\n");}
