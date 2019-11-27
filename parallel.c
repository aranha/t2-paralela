#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

/* CONTANTES */
#define GRAU 400
#define TAM_INI 1000000
#define TAM_INC 1000000
#define TAM_MAX 10000000

/* VARIAVEIS GLOBAIS */
double x[TAM_MAX], y[TAM_MAX], gabarito[TAM_MAX];

/* PROTOTIPOS */
double polinomio(double v[], int grau, double x);
void erro(char *msg_erro);

double polinomio(double a[], int grau, double x) {
  int i;
  double res = a[0], pot = x;
  for (i=1; i<=grau; ++i) {
    res+=a[i]*pot;
    pot=pot*x;
  }
  return res;
}

void erro(char *msg_erro) {
  fprintf(stderr, "ERRO: %s\n", msg_erro);
  MPI_Finalize();
  exit(1);
}

int main(int argc, char **argv) {
  int id; /* Identificador do processo */
  int n;  /* Numero de processos */
  int i, size;
  double *vet, valor, *vresp, resposta, tempo, a[GRAU + 1];
  int hostsize; /* Tamanho do nome do nodo */
  char hostname[MPI_MAX_PROCESSOR_NAME];
  MPI_Status status; /* Status de retorno */

  int tamanhoBloco; 
  int inicioBloco; 

  MPI_Init(&argc, &argv);
  MPI_Get_processor_name(hostname, &hostsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &n);
  
  if (id == 0) {
    /* Gera os coeficientes do polinomio */
    #pragma omp parallel for
    for (i=0; i<=GRAU; ++i)
      a[i] = (i%3==0)?-1.0:1.0;
    //distribui os coeficientes por broadcast para os escravos
    MPI_Bcast(a, GRAU + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    /* Preenche vetores */
    #pragma omp parallel for
    for (i=0; i<TAM_MAX; ++i) {
      x[i]=0.1+0.1*(double)i/TAM_MAX;
      gabarito[i]=polinomio(a, GRAU, x[i]);
    }
  } else {
    //escravo recebe os coeficientes do mestre
    MPI_Bcast(a, GRAU + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  /* Gera tabela com tamanhos e tempos */
  for (size=TAM_INI; size<=TAM_MAX; size+=TAM_INC) {
    //divide o tamanho do bloco de forma igual aos escravos
    tamanhoBloco = size/(n-1);
    inicioBloco = 0;

    if (id==0) {
        tempo = -MPI_Wtime();
        //distribuição de trabalho
        for (i = 1; i < n; ++i) {
          //envia o trabalho para o escravo (sub-array)
            MPI_Send(&x[inicioBloco], tamanhoBloco, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
            inicioBloco+=tamanhoBloco;
        }
        /* Calcula */
        for (i=inicioBloco; i<size; ++i){
            y[i]=polinomio(a, GRAU, x[i]);
        }
        inicioBloco=0;
        //recebe os resultados dos escravos
        for (i=1; i<n; ++i) {
            MPI_Recv(&y[inicioBloco], tamanhoBloco, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
            inicioBloco+=tamanhoBloco;
        }
        tempo+=MPI_Wtime();
        /* Verificacao */
        for (i = 0; i < size; ++i) {
            if (y[i]!=gabarito[i]) {
                erro("verificacao falhou!");
            }
        }
        /* Mostra tempo */
        printf("%d %lf\n", size, tempo);
    } else {
        //escravo recebe o trabalho, calcula e devolve o resultado ao mestre
        MPI_Recv(x, tamanhoBloco, MPI_DOUBLE, 0, id, MPI_COMM_WORLD, &status);
        for (i=0; i<tamanhoBloco; ++i){
            y[i]=polinomio(a, GRAU, x[i]);
        }
        MPI_Send(y, tamanhoBloco, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
    }
  }
  MPI_Finalize();
  return 0;
}