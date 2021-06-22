#include<stdio.h>
#include<math.h>
#include<time.h>

#define generation	1000
#define PC		0.25
#define PM		0.01

#define NPOP		10
#define BITCHROM 	39
#define m1		21
#define m2		18
#define PI		3.142857
#define TRUE		1
#define FALSE		0
#define WAIT	   	1	

/* Return a random integer from low to high */

short rrandom (low, high)
short low, high;
{
	
	long random();
        short nb = high - low + 1;
        return(random() % nb + low);
}

/* Return a random real from (0.0 1.0).*/

double ssrand ()
{
   short rrandom();

   return((double)(1 + rrandom(1, 31998)) / 32000.0);
}

/* Return a random bit (0 or 1).*/

int brand ()
{
   short rrandom();

   return(rrandom(0, 1));
}

int main(int argc, char **argv)
{
	int i,j,k,n,strong,findparent,cutpoint,regeneration,generationresult;
	int gentemp[1][BITCHROM],gene[NPOP][BITCHROM],gensl[NPOP][BITCHROM],genco[NPOP][BITCHROM],genmutt[NPOP][BITCHROM];
	int selected[NPOP],decx1[NPOP],decx2[NPOP],parent[1];
	float rk,totfx,fxresult,x1[NPOP],x2[NPOP],fx[NPOP],pk[NPOP],q[NPOP],r[NPOP];

	fxresult=0;	
	regeneration=0;
	
	for(;;)
	{
	
	time_t seconds;
	time(&seconds);
	srand((unsigned int) seconds);

	totfx=0;	
	strong=0;
        
	printf("\nGeneration - %d\n\n",regeneration);

	//Initial Population
	for (i=0;i<NPOP;i++)	
	{
		printf("v%d=",i);
		for (j=1;j<=BITCHROM;j++) 
		{
			if (regeneration<1)
			{
				gene[i][j]=brand();
			} else {
				gene[i][j]=genmutt[i][j];
			}
			printf("%d",gene[i][j]);
		}
		printf("\n");
	}

	printf("\n");
	
	//Evaluation
	for (i=0;i<NPOP;i++)
	{
	
		decx1[i]=0;
		n=m1-1;
		for (j=1;j<=m1;j++,n--)
		{
			decx1[i]=decx1[i]+gene[i][j]*pow(2,n);
			//printf("%d",gene[i][j]);
		}
		x1[i]=-3.0+decx1[i]*((12.1-(-3))/(pow(2,m1)-1));
		//printf("=%d, x1=%f\n",decx1[i],x1[i]);
		
		decx2[i]=0;
		n=m2-1;
		for (j=m1+1;j<=m1+m2;j++,n--)
		{
			decx2[i]=decx2[i]+gene[i][j]*pow(2,n);
			//printf("%d",gene[i][j]);
		}
		x2[i]=-3.0+decx2[i]*((12.1-(-3))/(pow(2,m2)-1));
		//printf("=%d, x2=%f\n",decx2[i],x2[i]);
		
		fx[i]=21.5 + x1[i] * sin(4*PI*x1[i]) + x2[i] * sin(20*PI*x2[i]);
		printf("fx%d(%f,%f)=%f\n",i,x1[i],x2[i],fx[i]);
	}
	
	for (i=0;i<NPOP;i++)
	{
		if (fx[strong]<fx[i])
		{
			strong=i;
		}
	}
	printf("strongest fx%d=%f\n",strong,fx[strong]);

	if (fxresult<=fx[strong]) 
	{
		fxresult=fx[strong];
		generationresult=regeneration;
	}
	
	//Selection
	for (i=0;i<NPOP;i++)
	{
		totfx=totfx+fx[i];
	}
	printf("totfx=%f\n",totfx);

	printf("\n");

	q[-1]=0;
	for (i=0;i<NPOP;i++)
	{
		pk[i]=fx[i]/totfx;
		q[i]=q[i-1]+pk[i];
		r[i]=ssrand();
		printf("pk(%d)=%f      q(%d)=%f r(%d)=%f\n",i,pk[i],i,q[i],i,r[i]);
	}

	for (i=0;i<NPOP;i++)
	{
		for (j=1;j<NPOP;j++)
		{
			if (q[j]>=r[i])
			{
				selected[i]=j;
				break;
			}
		} 
		printf("select[%d]=%d\n",i,selected[i]);
	}
	for (i=0;i<NPOP;i++)
	{
		printf("vsel%d from %d=",i,selected[i]);
		for (j=1;j<=BITCHROM;j++)
		{
			gensl[i][j]=gene[selected[i]][j];
			printf("%d",gensl[i][j]);
		}
		printf("\n");
	}

	printf("\n");
	
	//Crossover
	findparent=TRUE;
	j=0;i=0;
	while (findparent) 
	{
		r[i]=ssrand();
		if (r[i]<PC)
		{
			parent[j]=i;
			j++;
		}

		if (i>NPOP-1) { 
		 	i=0; } 
		else if (i<NPOP-1) {
		 	i++; }
		if (j>2) findparent=FALSE;
	}

	printf("father %d\n",parent[0]);
	printf("mother %d\n",parent[1]);

	for (i=0;i<NPOP;i++)
	{
		
		for (j=1;j<=BITCHROM;j++)
		{
			genco[i][j]=gensl[i][j];
		}
	}
	
	cutpoint=rrandom(1,BITCHROM);
	printf("cutpoint=%d\n",cutpoint);
	for (i=cutpoint;i<=BITCHROM;i++)
	{
		gentemp[0][i]=genco[parent[0]][i];
		genco[parent[0]][i]=genco[parent[1]][i];
		genco[parent[1]][i]=gentemp[0][i];
	}

	for (i=0;i<NPOP;i++)
	{
		printf("vco%d=",i);
		for (k=1;k<=BITCHROM;k++)
		{
			printf("%d",genco[i][k]);
		}
		printf("\n");
	}

	printf("\n");
	
	//Mutation
	for (i=0; i<NPOP; i++)
	{
		for (j=1;j<=BITCHROM;j++)
		{
			genmutt[i][j]=genco[i][j];
		}
		
	}

	for (i=0;i<NPOP;i++)
	{
		for (j=1;j<=BITCHROM;j++)
		{
			rk=ssrand();
			//printf("%f\n",rk);
			if (rk<PM) 
			{
				switch (genmutt[i][j]) {
					case 0: genmutt[i][j]=1;
						printf("gene 0->1 ");
						break;
					case 1: genmutt[i][j]=0;
						printf("gene 1->0 ");
						break;
				}
				
				printf("gene mutated [%d][%d] by random %f\n",i,j,rk); 		
			}
		}
	}

	for (i=0;i<NPOP;i++)
	{
		printf("vmu%d=",i);
		for (j=1;j<=BITCHROM;j++)
		{
			printf("%d",genmutt[i][j]);
		}
		printf("\n");
	}

	printf("\nGeneration %d, the strongest fx%d=%f\n",regeneration,strong,fx[strong]);
	printf("Result in %d generation with fx=%f\n",generationresult,fxresult);
	regeneration++;
	if (regeneration>generation) break;
	sleep(WAIT);
	}
}
