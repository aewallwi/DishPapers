#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define Nsample 121
main()
{
int i;
float freq[Nsample];
float feed_DB[Nsample],feed_P[Nsample];
float HERA_DB[Nsample],HERA_P[Nsample];
FILE *fp1=NULL;
FILE *fp2=NULL; 
FILE *fp3=NULL;
FILE *fp4=NULL;
FILE *fp5=NULL;
FILE *fp6=NULL;

fp1= fopen("useHERA_feeds11.csv","r");
if(fp1==NULL)
printf("Can't open the file");
fp2= fopen("simulation_feeds11_DB.csv","w");
if(fp2==NULL)
printf("Can't open the file");
fp3= fopen("simulation_feeds11_P.csv","w");
if(fp3==NULL)
printf("Can't open the file");


fp4= fopen("useHERAs11.csv","r");
if(fp4==NULL)
printf("Can't open the file");
fp5= fopen("simulation_HERAs11_DB.csv","w");
if(fp5==NULL)
printf("Can't open the file");
fp6= fopen("simulation_HERAs11_P.csv","w");
if(fp6==NULL)
printf("Can't open the file");



for(i=0;i<Nsample;i++)
{
fscanf(fp1,"%f\t\t%f\t\t%f\n",&freq[i],&feed_DB[i],&feed_P[i]);
fscanf(fp4,"%f\t\t%f\t\t%f\n",&freq[i],&HERA_DB[i],&HERA_P[i]);
fprintf(fp2,"%f,%10.9f\n", freq[i],feed_DB[i]);
fprintf(fp3,"%f,%10.9f\n", freq[i],feed_P[i]);
fprintf(fp5,"%f,%10.9f\n", freq[i],HERA_DB[i]);
fprintf(fp6,"%f,%10.9f\n", freq[i],HERA_P[i]);

//printf("Its running it new\n");
//printf("%d,%f\n", i,0.0);
}
fclose(fp1);
fclose(fp2);
//fclose(fp3);
//fclose(fp4);
//fclose(fp5);
//fclose(fp6);
}