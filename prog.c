#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int selection = 0;
float E = 1;
float P0 = 2;
float P1 = 1;

int remplirTab(char *fic, float **tableau, int n)
{
    FILE *tabfile = fopen(fic, "w+");
    for (int i = 0; i <= n; i++)
    {
        fprintf(tabfile, "%f %f %f\n", tableau[i][0], tableau[i][1], tableau[i][2]);
    }
    fclose(tabfile);
}

float f(int i, float yn, float zn, int tn)
{

    if (i == 1)
    {
        return zn;
    }
    else
    {
        if (selection == 2)
        {
            return (-P0 * tn * tn + P1 - E) * yn;
        }
        else
        {
            return yn * (-E);
        }
    }
}

float An(int i, float h, float yn, float zn, float tn)
{
    return h * f(i, yn, zn, tn);
}

float Bn(int i, float h, float yn, float zn, float tn, float A1, float A2)
{
    return h * f(i, yn + A1 / 2, zn + A2 / 2, tn + h / 2);
}

float **rk2(int N, float h, float yn, float zn, float tn, float **graphRK2)
{

    float A1 = 0;
    float A2 = 0;
    float B1 = 0;
    float B2 = 0;

    for (int i = 1; i <= N; i++)
    {
        A1 = An(1, h, yn, zn, tn);
        A2 = An(2, h, yn, zn, tn);
        B1 = Bn(1, h, yn, zn, tn, A1, A2);
        B2 = Bn(2, h, yn, zn, tn, A1, A2);

        float yn1 = yn + B1;
        float zn1 = zn + B2;

        yn = yn1;
        zn = zn1;
        tn += h;

        // printf("Yn : %f, Zn : %f, x: %f \n",  yn1, zn1, tn);

        graphRK2[i][0] = tn;
        graphRK2[i][1] = yn1;
        graphRK2[i][2] = zn1;
    }

    return graphRK2;
}

float **rk4(int N, float h, float yn, float zn, float tn, float **graphRK4)
{

    float A1 = 0;
    float A2 = 0;
    float B1 = 0;
    float B2 = 0;
    float C1 = 0;
    float C2 = 0;
    float D1 = 0;
    float D2 = 0;

    for (int i = 1; i <= N; i++)
    {
        A1 = An(1, h, yn, zn, tn);
        A2 = An(2, h, yn, zn, tn);
        B1 = Bn(1, h, yn, zn, tn, A1, A2);
        B2 = Bn(2, h, yn, zn, tn, A1, A2);
        C1 = Bn(1, h, yn, zn, tn, B1, B2);
        C2 = Bn(2, h, yn, zn, tn, B1, B2);
        D1 = Bn(1, h, yn, zn, tn, C1, C2);
        D2 = Bn(2, h, yn, zn, tn, C1, C2);

        float yn1 = yn + (A1 + 2 * B1 + 2 * C1 + D1) / 6;
        float zn1 = zn + (A2 + 2 * B2 + 2 * C2 + D2) / 6;

        yn = yn1;
        zn = zn1;
        tn += h;

        // printf("Yn : %f, Zn : %f, x: %f \n",  yn1, zn1, tn);

        graphRK4[i][0] = tn;
        graphRK4[i][1] = yn1;
        graphRK4[i][2] = zn1;
    }

    return graphRK4;
}

float **Euler(int N, float h, float y, float z, float a, float **tableau)
{
    float y1, y2;
    float hn = 0;

    for (int i = 1; i <= N; i++)
    {
        tableau[i][0] = tableau[i - 1][0] + h;
        tableau[i][1] = tableau[i - 1][1] + h * tableau[i - 1][2];
        tableau[i][2] = tableau[i - 1][2] - h * (tableau[i - 1][2] + hn * tableau[i - 1][1]);
        hn += h;
    }
    return tableau;
}

int main()
{
    int n;
    int correct = 1;
    int choixMethode;
    FILE *gnupipe = NULL;
    char *CommandeGnu[] = {"set terminal jpeg size 800,600", "set output 'Courbe.jpeg'", "unset key", "set title \"Courbe\"", "plot 'coordonnee.txt' with linespoints"};

    printf("Bonjour, pour résoudre les équations obtenues pendant l'étape de la modélisation, veuillez choisir l'une des méthodes suivant : \n1-Méthode RK2\n2-Méthode RK4\n3-Méthode d'EULER\nVotre choix : ");
    scanf("%d", &choixMethode);
    puts("\n");

    while (choixMethode != 1 && choixMethode != 2 && choixMethode != 3)
    {
        printf("Il y a une erreur dans la saisie de la méthode, veuillez entrer une réponse correcte.\n");
        scanf("%d", &choixMethode);
    }

    printf("Veuillez rentrer le nombre de points que vous souhaitez calculer : ");
    scanf("%d", &n);

    while (selection != 1 && selection != 2 && selection != 3)
    {
        if (!correct)
        {
            printf("Erreur dans la saisie, veuillez recommencer");
        }
        printf("Veuillez choisir un intervalle pour le potentiel: \n 1 pour x∈ [-3,-1] \n 2 pour x∈ [-1,1] \n 3 pour x∈ [1,3]\n => ");
        scanf("%d", &selection);
        if (selection != 1 && selection != 2 && selection != 3)
        {
            correct = 0;
        }
    }

    printf("Veuillez entrer une valeur pour E: ");
    scanf("%f", &E);

    P0 = 2 * E;
    P1 = 3 * E;

    int a, b;
    switch (selection)
    {
    case 1:
        a = -3;
        b = -1;
        break;
    case 2:
        a = -1;
        b = 1;
        break;
    case 3:
        a = 1;
        b = 3;
        break;
    default:
        a = -1;
        b = 1;
        printf("défaut");
    }

    float y0;
    printf("Veuillez rentrer la valeur de la fonction en %d : ", a);
    scanf("%f", &y0);
    float z0;
    printf("Veuillez entrer la valeur de la dérivé de la fonction en %d : ", a);
    scanf("%f", &z0);

    float h = ((float)b - (float)a) / ((float)n - 1);

    float **tableau = (float **)malloc((n + 1) * sizeof(float *));
    for (int i = 0; i <= n; i++)
    {
        tableau[i] = (float *)malloc(3 * sizeof(float));
    }

    tableau[0][0] = a;
    tableau[0][1] = y0;
    tableau[0][2] = z0;

    switch (choixMethode)
    {
    case 1:
        tableau = rk2(n, h, y0, z0, a, tableau);
        break;

    case 2:
        tableau = rk4(n, h, y0, z0, a, tableau);
        break;

    case 3:
        tableau = Euler(n, h, y0, z0, a, tableau);
        break;

    default:
        printf("erreur dans la résolution de l'équation... \n");
        break;
    }

    remplirTab("coordonnee.txt", tableau, n);
    gnupipe = popen("gnuplot", "w");
    for (int i = 0; i < 5; i++)
    {
        fprintf(gnupipe, "%s\n", CommandeGnu[i]);
    }
}
