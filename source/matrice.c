#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrice.h"


void annulle_matrice(double *tab, int lignes, int colonnes)
{
  int i, j;

  for (i=0; i < lignes; i++)
  {
    for (j=0; j < colonnes; j++)
    {

      *(tab + i * colonnes + j) = 0.0;

    } /* end for (j=0; j < colonnes; j++) */

  } /* end for (i=0; i < lignes; i++) */

} /* end void annulle_matrice(double *tab, int lignes, int colonnes) */

 /* D�termine la matrice inverse d'une matrice carr�e non nulle. La m�thode
    exploit�e est la d�composition en deux matrices triangulaires (L - R)
	qui sont chacune facilement inversable et dont le produit donne la 
	matrice inverse cherch�e.

    Cette m�thode souffre du d�faut qu'il faut que les mineurs de la diagonale
	ne soient pas null. Pour palier � cet inconv�nient, on introduit au d�part
	une l�g�re modification de la matrice (ajout d'une petite valeur � tous les
	�l�ments de la diagonale, ce qui permet de toujours calculer. 

    Cette erreur volontaire est minimis�e ensuite en corrigeant la matrice obtenue
	par it�rations jusqu'a obtenir une erreur inf�rieure � un miliardi�me
	de l'erreur introduite au d�part.

    Cette fa�on de faire n'est pas tr�s scientifique, mais c'�tait �cris pour le
	fun et de plus cela marche.
   
    Les entr�es sont:
     - un pointeur sur le premier �l�ments de la matrice;
     - un entier donnant le nombre de ligne et de colonne de la matrice
       carr�e;
     - un pointeur sur le permier �l�ment de la matrice inverse. Cette 
       matrice doit �tre d�finie dans le code appelant la fonction.

    La fonction retourne le d�terminant approch� de la matrice � inverser.
	Si ce d�terminant est nul, la matrice inverse d�termin�e est fausse (ou
	probl�me d'allocation dynamique de m�moire)
 
 */

double inverse_matrice(double *mat, int dim, double *mat_inv)
{
  
  double *U; /* matrice de travail */
  double *L; /* matrice de travail */
  double *Y; /* matrice de travail */
  double *Z; /* matrice de travail */
  double somme, erreur_acceptee, test;
  int i,j,k,p;

  if((U = (double *) malloc(sizeof(double) * (dim*dim)))==NULL) {return(0.0);}
  if((L = (double *) malloc(sizeof(double) * (dim*dim)))==NULL) {return(0.0);}
  if((Y = (double *) malloc(sizeof(double) * (dim*dim)))==NULL) {return(0.0);}
  if((Z = (double *) malloc(sizeof(double) * (dim*dim)))==NULL) {return(0.0);}

  /* annullation des matrices de travail */
  
  annulle_matrice((double*)U, dim, dim);  
  annulle_matrice((double*)L, dim, dim);
  annulle_matrice((double*)mat_inv, dim, dim);
  annulle_matrice((double*)Y, dim, dim);
  annulle_matrice((double*)Z, dim, dim);
  
  /* 
     d�termination de l'ordre de grandeur de la matrice pour �viter
     les erreurs de division par z�ro. Lorsque sur la diagonale, se 
     trouve des z�ros ou si un m�me nombre forme un coin dans la matrice,
     il y a risque de division par z�ro. Pour annuller ce risque, la  
     fonction ajoute aux �l�ments de la diagonale un nombre �gal au
     mili�mme de la valeur absolue non nulle la plus base. Cette 
     proc�dure introduit une erreur de calcul de l'ordre du mili�me 
     mais permet de calculer la matrice inverse de n'importe qu'elle 
     matrice carr�e non nulle.
	 
	 le d�terminant est une valeur approch�e (au milli�me pr�s)
	 
	 la matrice est ensuite corrig�e par un algorhitme r�cursif
	 r�alis� au maximum 10 fois ou jusqu'� ce que l'erreur entre
	 la matrice unitaire et la multiplication de la matrice par son
	 inverse soit plus petit qu'un miliardi�me de l'erreur accept�e.
  */
  
  erreur_acceptee = 1000.0;
  somme = 0.0;
  for(k=0;k<dim;k++)      /* d�termination de la valeurs absolue non nulle */
  {                       /* la plus basse */
    for(i=0;i<dim;i++)
    {
    test = fabs(*(mat + k*dim +i));
    if(somme == 0.0)
      { 
        somme = test;
      } /* end if */
    if(test>0 && test< somme)
      {
        somme = test;
      } /* end if */
    } /* end for(i=0;i<dim;i++) */
  } /* end for(k=0;k<dim;k++) */
  erreur_acceptee = somme/erreur_acceptee;

  for(k=0;k<dim;k++)  /* modifie la diagonale de la matrice d'origine */
  {
    (*(mat+k*dim +k)) += erreur_acceptee; 
  }

  /* d�termination des param�tres de L et U */

  for (k=0; k<dim; k++)
  {
    for (j=k; j<dim; j++)
    {
      somme = 0.0;
      
      for(p=0; p<k; p++)
      {
  		  somme += (*(L+k*dim+p)) * (*(U+p*dim+j));
      } /* end for(p=0; p<k; p++) */

  	  *(U+k*dim+j) = (*(mat + k*dim + j)) - somme;
   
      somme = 0.0;
      
      for(p=0; p<k; p++)
      {
  		  somme += (*(L+j*dim+p)) * (*(U+p*dim+k));
      } /* end for(p=0; p<k; p++) */

  	  *(L+j*dim+k) = 1.0;
      
      if(j!=k)
      {
  		  *(L+j*dim+k) = ((*(mat + j*dim + k)) - somme)/(*(U+k*dim+k));
      }

    } /* end for (j=k; j<dim; j++) */
  } /* end for (k=0; k<dim; k++) */

  for(k=0;k<dim;k++)  /* restitue la matrice d'origine */
  {
    (*(mat+k*dim+k)) -= erreur_acceptee;  
  }

  /* inverse les matrices L et U */

  for (j=0; j<dim; j++)
  {
  	  *(Y+j*dim+j)=1.0;

    for (i=j+1; i<dim; i++)
    {
      somme = 0.0;
      
      for(k=j; k<i; k++)
      {
  		  somme += (*(L+i*dim+k)) * (*(Y+k*dim+j));
      } /* end for(k=j; k<i; k++) */

  	  *(Y+i*dim+j) = (0.0 - somme)/(*(L+i*dim+i));

    } /* end for (i=j+1; i<dim; i++) */
  } /* end for (j=0; j<dim; j++) */

  for (j=0; j<dim; j++)
  {
  	  *(Z+j*dim+j)=1.0/(*(U+j*dim+j));

    for (i=j-1; i>-1; i--)
    {
      somme = 0.0;
      
      for(k=i+1; k<j+1; k++)
      {
  		somme += (*(U+i*dim+k)) * (*(Z+k*dim+j));
      } /* end for(k=j; k<i-1; k++) */

  	  *(Z+i*dim+j) = (0.0 - somme)/(*(U+i*dim+i));

    } /* end for (i=j+1; i<dim; i++) */
  } /* end for (j=0; j<dim; j++) */
  

/* calcul de la matrice inverse par multiplication de Z par Y */

  for (i=0; i<dim; i++)
  {  
    for (j=0; j<dim; j++)
      {
        for (k=0; k<dim; k++)  
        {
  			(*(mat_inv + i*dim +j)) += (*(Z+i*dim+k)) * (*(Y+k*dim+j));
        } /* end for (k=0; k<dim; k++) */
      } /* end for (j=0; j<dim; j++) */
  } /* end for (i=0; i<dim; i++) */

  /* 
     calcul du d�terminant de la matrice par multiplication des
     �l�ments de la diagonale de U 
  */

  somme = 1.0;
  for(k=0;k<dim;k++)
  {
  	somme *= (*(U+k*dim+k));
  } /* end for(k=0;k<dim;k++) */
  
  /* l'introduction de l'erreur_accept�e ne permet plus de trouver un
     d�terminant nul. Le test suivant impose un d�terminant nul
     si celui-ci est de l'ordre de grandeur de l'erreur_accept�e */
  
  if(fabs(somme) < erreur_acceptee*10) 
    {
      somme = 0.0;
    } /* end if(somme < erreur_acceptee*10) */
	
/* correction de la matrice inverse si le d�terminant n'est pas nul */

if (somme != 0.0)
  {

  annulle_matrice((double*)U, dim, dim);  
  annulle_matrice((double*)L, dim, dim);
  annulle_matrice((double*)Y, dim, dim);
  annulle_matrice((double*)Z, dim, dim);

  /* creation de la matrice unite U */
  for(k=0;k<dim;k++)
    {
  	  *(U+k*dim+k) = 1.0;
    } /* end for(k=0;k<dim;k++) */

    /* L est le resultat de la multiplication de la matrice par son inverse approchee */
    mutiplication_matrice((double*)mat, dim, dim, dim, (double*)mat_inv, (double*)L);
   /* Y est le resultat de la matrice unite (U) moins la matrice L */
    soustraction_matrice((double*)U, dim, dim, (double*)L, (double*)Y);
  
  k = 0;    
  	  test = 0.0;
      for(i=0;i<dim;i++)
       {
	     for(j=0;j<dim;j++)
		  {
			if(test<fabs(*(Y+i*dim+j))) {test = fabs(*(Y+i*dim+j));}
		  } /* end for(j=0;j<dim;j++) */
       } /* end for(i=0;i<dim;i++) */
	   
	   printf("pour k = %d, test = %e\n", k, test);

  do
    {
  
       /* Z est le resultat de la multiplication de la matrice inverse par Y */
      mutiplication_matrice((double*)mat_inv, dim, dim, dim, (double*)Y, (double*)Z);
  
      /* la matrice inverse corrig�e est la somme de la matrice inverse et de Z */
      addition_matrice((double*)mat_inv, dim, dim, (double*)Z, (double*)mat_inv);
	  
	  /* test pour la fin du traitement sur base de la difference entre U et mat * mat_inv 
	     et preparation du calcul pour le tour suivant si necessaire */
    
	  /* L est le resultat de la multiplication de la matrice par son inverse approchee */
      mutiplication_matrice((double*)mat, dim, dim, dim, (double*)mat_inv, (double*)L);
  	  
     /* Y est le resultat de la matrice unite (U) moins la matrice L */
      soustraction_matrice((double*)U, dim, dim, (double*)L, (double*)Y);
	  
	  test = 0.0;
      for(i=0;i<dim;i++)
       {
	     for(j=0;j<dim;j++)
		  {
			if(test<fabs(*(Y+i*dim+j))) {test = fabs(*(Y+i*dim+j));}
		  } /* end for(j=0;j<dim;j++) */
       } /* end for(i=0;i<dim;i++) */
	   
	   k++;
 	   
    } while (test>erreur_acceptee * 1e-9 && k<10);
	
  } /* end if (somme != 0.0) */

  free(U);
  free(L);
  free(Y);
  free(Z);

  return(somme); /* si la somme est nulle, la matrice inverse est fausse */

} /* end double inverse_matrice(double *mat, int dim, double *mat_inv) */

/*  
  Multiplication de deux matrices mat1*mat2 de dimentions:
    
    - mat1[ligne1][colonne1]
    - mat2[colonne1][colonne2]

  dans une matrice resultat[ligne1][colonne2]

  Les donn�es sont:

    - pointeur sur le premier �l�ment de la matrice 1
    - entier du nombre de lignes de la matrice 1
    - entier du nombre de colonnes de la matrice 1 = nombre de lignes de la
      matrice 2
    - entier du nombre de colonnes de la matrice 2
    - pointeur sur le premier �l�ment de la matrice 2
    - pointeur sur le premier �l�ment de la matrice de resultat

*/

void mutiplication_matrice(double *mat1, int ligne1, int colonne1, int colonne2, double *mat2, double *resultat)
{
  int i, j, k;  

  annulle_matrice((double*)resultat, ligne1, colonne2); 

  for (i=0; i<ligne1; i++)
  {  
    for (j=0; j<colonne2; j++)
      {
        for (k=0; k<colonne1; k++)  
        {
          /* (*(resultat + i*ligne1 +j)) += (*(mat1 + i*ligne1 + k)) * (*(mat2 + k*colonne1 +j)); !!! faux faux pour matrice non carr�e */
			(*(resultat + i*colonne2 +j)) += (*(mat1 + i*colonne1 + k)) * (*(mat2 + k*colonne2 +j)); /* RR2.0 : 21/06/2012 : correction bug */
        } /* end for (k=0; k<dim; k++) */
      } /* end for (j=0; j<dim; j++) */
  } /* end for (i=0; i<dim; i++) */


} /* end void multiplication_matrice(double *mat, int dim, double *mat_inv, double *unitaire) */
 
/*  
  addition de deux matrices mat1 et mat2 de dimentions:
    
    - mat1[lignes][colonnes]
    - mat2[lignes][colonnes]

  dans une matrice resultat[lignes][colonnes]

  Les donn�es sont:

    - pointeur sur le premier �l�ment de la matrice 1
    - entier du nombre de lignes des matrices
    - entier du nombre de colonnes des matrices
    - pointeur sur le premier �l�ment de la matrice 2
    - pointeur sur le premier �l�ment de la matrice de resultat

*/

void addition_matrice(double *mat1, int lignes, int colonnes, double *mat2, double *resultat)
{

int i, j;

  for (i=0; i<lignes; i++)
  {  
    for (j=0; j<colonnes; j++)
      {
          /* (*(resultat + i*lignes +j)) = (*(mat1 + i*lignes + j)) + (*(mat2 + i*lignes +j)); !!!!! faux pour matrice non carr�e */
		  (*(resultat + i*colonnes +j)) = (*(mat1 + i*colonnes + j)) + (*(mat2 + i*colonnes +j)); /* RR2.0 : 21/06/2012 : correction bug */
      } /* end for (j=0; j<dim; j++) */
  } /* end for (i=0; i<dim; i++) */



} /* end void addition_matrice(double *mat1, int lignes, int colonnes, double *mat2, double *resultat) */
 
/*  
  soustraction de deux matrices mat1 - mat2 de dimentions:
    
    - mat1[lignes][colonnes]
    - mat2[lignes][colonnes]

  dans une matrice resultat[lignes][colonnes]

  Les donn�es sont:

    - pointeur sur le premier �l�ment de la matrice 1
    - entier du nombre de lignes des matrices
    - entier du nombre de colonnes des matrices
    - pointeur sur le premier �l�ment de la matrice 2
    - pointeur sur le premier �l�ment de la matrice de resultat

*/
 
void soustraction_matrice(double *mat1, int lignes, int colonnes, double *mat2, double *resultat)
{

int i, j;

  for (i=0; i<lignes; i++)
  {  
    for (j=0; j<colonnes; j++)
      {
         /* (*(resultat + i*lignes +j)) = (*(mat1 + i*lignes + j)) - (*(mat2 + i*lignes +j)); !!!!! faux pour matrice non carr�e */
		  (*(resultat + i*colonnes +j)) = (*(mat1 + i*colonnes + j)) - (*(mat2 + i*colonnes +j)); /* RR2.0 : 21/06/2012 : correction bug */
      } /* end for (j=0; j<dim; j++) */
  } /* end for (i=0; i<dim; i++) */


} /* end void soustraction_matrice(double *mat1, int lignes, int colonnes, double *mat2, double *resultat)*/

/*  
  copie d'une matrice dans une autre deux matrices mat_copie - mat de dimentions identiques:
    
    - mat_copie[ligne][colonne]
    - mat2[ligne][colonne]

  
  Les donn�es sont:

    - pointeur sur le premier �l�ment de la matrice recevant la copie
	- pointeur sur le premier �l�ment de la matrice � copier
    - entier du nombre de lignes des matrices
    - entier du nombre de colonnes des matrices
    
*/


void copie_matrice(double *mat_copie, double *mat, int ligne, int colonne)
{
	int i, j;

	for (i=0; i<ligne; i++)
	{  
		for (j=0; j<colonne; j++)
		{
			/* (*(mat_copie + i*ligne +j)) = (*(mat + i*ligne + j)); !!!!! faux pour matrice non carr�e */
		  (*(mat_copie + i*colonne +j)) = (*(mat + i*colonne + j)); /* RR2.0 : 21/06/2012 : correction bug */
		} /* end for (j=0; j<dim; j++) */
	} /* end for (i=0; i<dim; i++) */
}

/*  
  multiplie une matrice par un r�el:
          
  Les donn�es sont:

    - pointeur sur le premier �l�ment de la matrice � multiplier
	- un double = r�el multiplicateur
    - entier du nombre de lignes de la matrice
    - entier du nombre de colonnes de la matrice
    
*/

void reelfoismatrice(double *mat, double reel, int ligne, int colonne) /* multiplie une matrice par un r�el, le r�sultat est dans la matrice de d�part */
{
	int i, j;

	for (i=0; i<ligne; i++)
	{  
		for (j=0; j<colonne; j++)
		{
			/* (*(mat + i*ligne +j)) *= reel; !!!!! faux pour matrice non carr�e */
		   (*(mat + i*colonne +j)) *= reel; /* RR2.0 : 21/06/2012 : correction bug */
		} /* end for (j=0; j<dim; j++) */
	} /* end for (i=0; i<dim; i++) */

}

/*  
  addition de deux matrices mat et mat_add de dimentions identiques:
    
    - mat[lignes][colonnes]
    - mat_add[lignes][colonnes]

  le r�sultat �tant plac� dans la matrice mat[lignes][colonnes]

  Les donn�es sont:

    - pointeur sur le premier �l�ment de la matrice 1
    - entier du nombre de lignes des matrices
    - entier du nombre de colonnes des matrices
    - pointeur sur le premier �l�ment de la matrice 2
    
*/


void addition_dans_matrice(double *mat, int lignes, int colonnes, double *mat_add) /* addition de mat+mat_add avec r�sultats dans mat */
{
  int i, j;

  for (i=0; i<lignes; i++)
  {  
    for (j=0; j<colonnes; j++)
      {
          /* (*(mat + i*lignes +j)) += (*(mat_add + i*lignes + j)); !!!!! faux pour matrice non carr�e */
		   (*(mat + i*colonnes +j)) += (*(mat_add + i*colonnes + j)); /* RR2.0 : 21/06/2012 : correction bug */
      } /* end for (j=0; j<dim; j++) */
  } /* end for (i=0; i<dim; i++) */



}

/* inversion de la matrice � l'aide des co�fficients d'un polynome caract�ristique:

  recherche du polynome (p) par la m�thode de Leverrier, 

          1
   p  = - - * (s  + p  * s    + .... + p    * s )
    n     n     n    1    n-1           n-1    1

	
  puis, calcul de la matrice inverse par la formule:

   -1    1     n-1         n-2
  A  = - - * (A    + p  * A     + ..... + p   * E)
         p            1                    n-1
          n

  Les donn�es sont:

     - un pointeur sur le premier �l�ments de la matrice;
     - un entier donnant le nombre de ligne et de colonne de la matrice
       carr�e;
     - un pointeur sur le permier �l�ment de la matrice inverse. Cette 
       matrice doit �tre d�finie dans le code appelant la fonction.

  retourne 1 si matrice est invers�e ou 0 si d�faut d'allocation dynamique de m�moire
  ou si un diviseur est nul

 */

int inverse_matrice_poly(double *mat, int dim, double *mat_inv)
{
	double *tr1, *tr2, *s, *p;
	int i, j;

	if((s = (double *) malloc(sizeof(double) * (dim+1)))==NULL) {return(0);} /* interm�diaires pour le calcul de p */
	if((p = (double *) malloc(sizeof(double) * (dim+1)))==NULL) {return(0);} /* co�fficients du polynome caract�ristique de mat */
	if((tr1 = (double *) malloc(sizeof(double) * (dim*dim)))==NULL) {return(0);} /* matrices de travail */
	if((tr2 = (double *) malloc(sizeof(double) * (dim*dim)))==NULL) {return(0);} /* matrices de travail */

	annulle_matrice((double*)tr1, dim, dim);
	annulle_matrice((double*)tr2, dim, dim);

	/* calcul des co�fficients du polynome caract�ristiques
	par le m�thode de Leverrier(avec les matrices de puissance -> dim */

	/* les co�ficients interm�diaires s sont les sommes des �l�ments des
	diagonales des matrices de puissance 1 � dim-1 */

	/* calcul de s1 */
	*(s) = 0.0;

	for (i=0; i < dim; i++)
	{
		 *(s) += *(mat + i * dim + i);
		  //printf("mat = %f\n", *(mat + i * dim + i));			  	    
	} /* end for (i=0; i < dim; i++) */

	copie_matrice(tr1, mat, dim, dim);

	/* calculs de s2 � s_dim */

	for(i=1; i<dim;i++)
	{
		mutiplication_matrice(tr1, dim, dim, dim, mat, tr2);
		*(s+i) = 0.0;

		for (j=0; j < dim; j++)
		{
			 *(s+i) += *(tr2 + j * dim + j);
		} /* end for (j=0; j < dim; j++) */

		copie_matrice(tr1, tr2, dim, dim);
	}

	/* d�termination des co�ficients p sur base des co�ficients s */

	*(p) = -(*(s));

	for(i=1; i<dim;i++)
	{
		*(p+i) = *(s+i);
		for(j = 0; j<i; j++)
		{
			*(p+i)+= (*(p+j)) * (*(s+i-j-1));
		}

		(*(p+i))/= -(i+1);
	}

	if(*(p+dim-1) == 0)
	{
		/* inversion matrice impossible */
		free(s);
		free(p);
		free(tr1);
		free(tr2);
		return(0);
	}

	annulle_matrice((double*)mat_inv, dim, dim);

	for(i=0; i<dim;i++) /* matrice unit� * p_n-1 */
	{
		*(mat_inv + i*dim + i) = *(p+dim-2);
	}

	copie_matrice(tr1, mat, dim, dim); /* copie mat dans tr1 */
	copie_matrice(tr2, mat, dim, dim); /* copie mat dans tr2 */

	for(i=2; i<dim; i++) /* calcul les �l�ments suivants */
	{
		reelfoismatrice(tr1, *(p+dim-(i+1)), dim, dim);
		addition_dans_matrice(mat_inv, dim, dim, tr1);
		copie_matrice(tr1, tr2, dim, dim); /* remet matrice puissance i dans tr1 avant mult par mat car tr1 modifier par fois reel */
		mutiplication_matrice(tr1, dim, dim, dim, mat, tr2);
		copie_matrice(tr1, tr2, dim, dim);		
	}

	addition_dans_matrice(mat_inv, dim, dim, tr1);
	reelfoismatrice(mat_inv, -1/(*(p+dim-1)), dim, dim);

	free(s);
	free(p);
	free(tr1);
	free(tr2);

	return(1); /* calcul termin� */

}

void print_matrice(double *tab, int lignes, int colonnes)
{
  int i, j;

  for (i=0; i < lignes; i++)
  {
    for (j=0; j < colonnes; j++)
    {

      printf("%f ", *(tab + i * colonnes + j));

    } /* end for (j=0; j < colonnes; j++) */

    printf("\n");

  } /* end for (i=0; i < lignes; i++) */

  printf("\n");

} /* end void print_matrice(double *tab, int lignes, int colonnes) */


/* inversion de la matrice par la m�thode du pivot.

  le pivot est choisi dans chaque ligne comme �tant le plus grand co�fficient non nul.

  La ligne de la matrice est divid�e par le pivot, puis cette ligne est retranch�e des
  autres lignes en la multipliant par le co�fficient de la ligne ayant le m�me indice
  que le pivot choisi. La colonne du pivot devient donc nulle sauf pour la ligne du pivot
  ou l'�l�ment est �gal � 1.

  On fait en parall�le la m�me chose dans le matrice unit� de d�part qui devient, en
  fin de calcul la matrice inverse dont les lignes sont m�lang�es. En effet, si pour la
  premi�re ligne de la matrice, c'est l'�l�ment x qui est le plus grand, le traitement
  donne le r�sultat pour la ligne x et non pas pour la premi�re ligne. Comme toutes les
  lignes vont �tre trait�e de 0 � n-1, les lignes de la matrice inverse seront m�lang�es.
  En gardant en m�moire l'ordre de traitement des �l�ments, il est facile en fin de
  calcul de la matrice inverse m�lang�es de recopier les lignes aux bon endroits dans
  la matrice inverse (dans notre cas, la premi�re ligne doit �tre copi�e dans la ligne
  x de la matrice inverse finale.

  Comme le syst�me d�truit la matrice d'entr�e, celle-ci est copi�e avant traitement.

  Le d�terminant est le produit de tous les pivots.

  Si une ligne ne contient plus que des �l�ments nuls, la matrice n'est pas inversable
  et son d�terminant est nul.
  
  Les donn�es sont:

     - un pointeur sur le premier �l�ments de la matrice;
     - un entier donnant le nombre de ligne et de colonne de la matrice
       carr�e;
     - un pointeur sur le permier �l�ment de la matrice inverse. Cette 
       matrice doit �tre d�finie dans le code appelant la fonction.

 */

double inverse_matrice_pivot(double *mat, int dim, double *mat_inv) /* d�termine la matrice inverse d'une matrice carr�e non nulle par la m�thode du pivot*/
{
	double *mtr, *melange; /* matrice de d�part de travail et matrice inverse m�lang�e */
	int *memoire;
	int i, j, k, indice;
	double determinant = 1.0;
	double pivot, diviseur;
	double max;

	if((mtr = (double *) malloc(sizeof(double) * (dim*dim)))==NULL) {return(0.0);} /* matrice de travail */
	if((melange = (double *) malloc(sizeof(double) * (dim*dim)))==NULL) {return(0.0);} /* matrice inverse m�lang�e */
	if((memoire = (int *) malloc(sizeof(int) * dim))==NULL) {return(0.0);}; /* tableau de correspondance entre la ligne et la colonne trait�e */
	copie_matrice(mtr, mat, dim, dim); /* copie mat dans mtr */
	annulle_matrice(melange, dim, dim); /* anulle tous les �l�ments de la matrice */
	annulle_matrice(mat_inv, dim, dim);
	
	for(i=0; i<dim;i++) /* matrice inverse m�lang�e commence comme matrice unit� */
	{
		*(melange + i*dim + i) = 1.0;
	}
	
	/* inversion de la matrice par recherche dans chaque ligne du pivot le plus grand en absolu */

	for(i=0;i<dim;i++) /* on traite la ligne i */
	{
		indice = 0;
		
		/* recherche de l'indice de l'�l�ment le plus grand en absolu 
		
		!!!! en prenant l'�l�ment le plus grand, on ne traite plus la ligne correspondant
		 � la i, mais la ligne correspondant � l'indice de l'�l�ment choisi. On construit donc la
		 matrice inverse avec des lignes m�lang�es. Il faut donc garder une trace de la 
		 correspondance entre la ligne trait�e et la ligne correspondante dans la matrice
		 inverse finale pour par la suite recopier le tout dans le bon ordre.
		
		*/
		
		max = fabs(*(mtr+i*dim+indice));
		j = 1;
		
		while(j<dim)
		{
			if(fabs(*(mtr+i*dim+j)) > max)
			{
				max = fabs(*(mtr+i*dim+j));
				indice = j;
			}		
			j++;
		}
		
		if(max == 0.0) /* a trouver une ligne dont tous les co�fficients sont nuls */
		{
			free(mtr);
			free(melange);
			free(memoire);
			return(0.0); /* matrice n'est pas inversible car d�terminant nul */
		}
		
		pivot = *(mtr+i*dim+indice);

		*(memoire+i) = indice; /* la ligne i correspond � la ligne finale indice */

		/* calcule le d�terminant */
		determinant*=pivot;
		
		/* divise la ligne des deux matrices par le pivot */
		for(j = 0; j<dim ; j++)
		{
			*(mtr+i*dim+j)/=pivot;
			*(melange+i*dim+j)/=pivot;
		}
		
		/* soustrait la ligne du pivot des autres lignes * le co�fficient de la colonne du pivot */

		for(j=0; j<i; j++)
		{
			diviseur = *(mtr+j*dim+indice);
		
			for(k=0; k<dim; k++)
			{
				*(mtr+j*dim+k)-= diviseur * (*(mtr+i*dim+k));
				*(melange+j*dim+k)-= diviseur * (*(melange+i*dim+k));
			}			
		}

		for(j=i+1; j<dim; j++)
		{
			diviseur = *(mtr+j*dim+indice);
			
			for(k=0; k<dim; k++)
			{
				*(mtr+j*dim+k)-= diviseur * (*(mtr+i*dim+k));
				*(melange+j*dim+k)-= diviseur * (*(melange+i*dim+k));
			}			
		}		
	}

	/* copie des lignes de la matrice inverse m�lang�e � la bonne place dans la matrice
	inverse finale : ligne i de m�lange = ligne *(memoire+i) de la matrice inverse  */

	for(i=0;i<dim;i++) /* on traite la ligne i */
	{
		for(j=0; j<dim; j++)
		{
			*(mat_inv+(*(memoire+i))*dim+j) = *(melange+i*dim+j);
		}

	}

	free(mtr);
	free(melange);
	free(memoire);

	return(determinant);
}
