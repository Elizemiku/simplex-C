/*Elizabeth Borgognoni Souto				RA: 170409
  Jessica Fernanda Teixeira do Prado			155875
  Larissa Couto Luiz							117557

Este programa implementa o algoritmo primal simplex,
distinguindo de acordo com a entrada já na forma padrão
qual a fase a ser implementada. 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/*funcoes utilizadas */

//funcao que escolhe as matrizes Basica, Nao basica, custosBasico e custoNbasico da fase 2 feita dentro da fase 1
void escolher_matrizBNF2(double **matrizA, double **matrizB, double **matrizN, int l, int c, double *vetCustos, double *vetCustoBasico, double *vetCustoNBasico, int escolha){
	
	int i, j, k;
	
	for(i = 0; i < l; i++){
		for(j = c-(c-l), k = 0; j < c && k < c - (c - l); j++, k++){
			matrizN[i][k] = matrizA[i][j];
		}
	}

	printf("\nMatriz N:\n");

	for(i = 0; i < l; i++){
		for(j = 0; j < c-l; j++){
			printf("%.2lf ", matrizN[i][j]);
		}
		printf("\n");
	}

	for(i = 0; i < l; i++){
		for(j = 0; j < c-(c-l); j++){
			matrizB[i][j] = matrizA[i][j];
		}
	}
	
	printf("\nMatriz B:\n");
	
	for(i = 0; i < l; i++){
		for(j = 0; j < c-(c-l); j++){
			printf("%.2lf ", matrizB[i][j]);
		}
		printf("\n");
	}

	for(i = c-(c-l), k = 0; i < c && k < c-l; i++, k++){
		vetCustoNBasico[k] = vetCustos[i];
	}
	
	printf("\n\nVetor custo nao basico:\n");

	for(i = 0; i < c-l; i++){
		printf("%.2lf ", vetCustoNBasico[i]);
	}

	for(i = 0; i < c-(c-l); i++){
		vetCustoBasico[i] = vetCustos[i];
	}

	printf("\nVetor custo basico:\n");

	for(i = 0; i < c-(c-l); i++){
		printf("%.2lf ", vetCustoBasico[i]);
	}

	printf("\n");
}

//funcao que escolhe as matrizes Basica, Nao basica, custosBasico e custoNbasico da fase 1
void escolher_matrizBNfase1(double **matrizA, double **matrizB, double **matrizN, int l, int c, double *vetCustos, double *vetCustoBasico, double *vetCustoNBasico){
	
	int i, j, k;
	
	for(i = 0; i < l; i++){
		for(j = c-l, k = 0; j < c && k < c - (c-l); j++, k++){
			matrizB[i][k] = matrizA[i][j];
		}
	}
	
	printf("\nMatriz B:\n");
	
	for(i = 0; i < l; i++){
		for(j = 0; j < c - (c-l); j++){
			printf("%.2lf ", matrizB[i][j]);
		}
		printf("\n");
	}

	for(i = 0; i < l; i++){
		for(j = 0; j < c-l; j++){
			matrizN[i][j] = matrizA[i][j];
		}
	}
	
	printf("\nMatriz N:\n");
	
	for(i = 0; i < l; i++){
		for(j = 0; j < c-l; j++){
			printf("%.2lf ", matrizN[i][j]);
		}
		printf("\n");
	}

	
	for(i = 0, k = c - l; i < c - l && k < c ; i++, k++){
		vetCustoBasico[i] = vetCustos[k];
	}

	printf("\nVetor custo basico:\n");

	for(i = 0; i < l; i++){
		printf("%.2lf ", vetCustoBasico[i]);
	}

	for(i = 0; i < c-l ; i++){
		vetCustoNBasico[i] = vetCustos[i];
	}


	printf("\n\nVetor custo nao basico:\n");	

	for(i = 0; i < c-l ; i++){
		printf("%.2lf ", vetCustoNBasico[i]);
	}
	printf("\n");

}

//funcao que realiza a troca os vetores basicos e nao basicos
void swapV(double *vetA, int swap1, int swap2){
	
	double aux;
	
	aux = vetA[swap1];
	vetA[swap1] = vetA[swap2];
	vetA[swap2] = aux;
}


//funcao que realiza a troca as matrizes basicas e nao basicas
void swapM(double **matrizA, int swap1, int swap2 , int l, int *trocas){
	
	double aux = 0;
	int i;
	
	for(i = 0; i < l; i++){
		aux = matrizA[i][swap1];
		matrizA[i][swap1] = matrizA[i][swap2];
		matrizA[i][swap2] = aux;
	}
	aux = trocas[swap1];
	trocas[swap1] = trocas[swap2];
	trocas[swap2] = aux;
}

//funcao que multiplica vetores 
double multiplica(double* vetL, double* colNB, int l){
	int i; 
	double som = 0;
	
	for(i = 0; i < l; i++){
		som = som + (vetL[i] * colNB[i]);
	}
	return som;
}

//funcao que calcula os custos relativos
int custos_relativos(double* vetCNB, double *vetLT, double **matrizN, int l, int c, int custozerado){
	
	int tamCNB = 0, i, j, indice = 0, cond = 1, cond1 = 0; 
	double *vetCol, menor = 0, *cr;
	
	tamCNB = c - l;
	vetCol = malloc(l * sizeof(double));
	cr = malloc((c-l) * sizeof(double));
	
	printf("\n");
	
	for(i = 0; i < tamCNB; i++){
		for(j = 0; j < l; j++){
			vetCol[j] = matrizN[j][i];
		}
		cr[i] = vetCNB[i] - multiplica(vetLT, vetCol, l);
		printf("CN%d = %.2lf \n", i + 1, cr[i]);
	}
	
	
	for(i = 0; i < tamCNB; i++){
		if(cr[i] < 0){
			cond = 0;
		}
		else if(cr[i] == 0){
			cond1++;
		}
	}
	
/*Solucoes multiplas sao possiveis quando um ou mais resultados do custo calculado e igual a 0 e alguns
dos meus vetores da solucao otima dao zero*/

/* A  existência  de  custos  relativos  nulos  e  condição  necessaria para multiplas solucoes otimas, mas 
nao e suficiente*/   

//analisa se alguns do custos calculados e igual a 0 
	if(cond == 1 && cond1 > 0){
		printf("\nPode ter multiplas solucoes otimas\n\n");
		custozerado = 1;
	}	

	if(cond == 0){
		for(i = 0; i < tamCNB; i++){
			if(cr[i] < menor){
				indice = i;
				menor = cr[i];
			}
		}
	
		printf("Menor indice = %d\n", indice + 1);
		printf("\nA coluna %d da matrizN entra na base!\n", indice + 1);

		free(vetCol);
		free(cr);
		return indice;
	}
	else if(cond == 1 && (cond1 != tamCNB)){
		free(vetCol);
		free(cr);
		return -1;
	}
	return 0;
}

//funcao que calcula a transposta de uma matriz
double** transpor_matriz(double **matriz, int l){
	double **transp; 
	int i, j;
	
	transp = malloc(l * sizeof(double));
	
	for(i = 0; i < l; i++){
		transp[i] = malloc(l * sizeof(double));
	}
	
	for(i = 0; i < l; i++){
		for(j = 0; j < l; j++){
			transp[j][i] = matriz[i][j];
		}
	}
	
	printf("\nMatriz Basica Transposta\n");
	
	for(i = 0; i < l; i++){
		for(j = 0; j < l; j++){
			printf("%lf ", transp[i][j]);
		}
		printf("\n");
	}
	
	return transp;
}

//funcao que libera memoria de uma matriz
void liberarMemoriaMatriz(double **matriz, int l){
	
	int i;
	
	for(i = 0; i < l; i++){
		free(matriz[i]);
	}
		
	free(matriz);
	
}


//funcao que resolve os sistemas calculados ao longo do algoritmo por meio de metodo de gauss de calculo numerico
void resolver_sistema(double **matrizB, int l, double *vetRecursos, double *vetSol){
	
    int k, i, j, posicao, erro;
	double m, soma, aux, pivo, **matrizCopia, *vetCopia;
	
    erro = 0; //quando da erro de divisao por zero

	matrizCopia = malloc(l * sizeof(double));
	
	for(i = 0; i < l; i++){
		matrizCopia[i] = malloc(l * sizeof(double));
	}
	
	for(i = 0; i < l; i++){
		for(j = 0; j < l; j++){
			matrizCopia[i][j] = matrizB[i][j];
		}
	}
	
	vetCopia = malloc(l * sizeof(double));
	
	for(i = 0; i < l; i++){
		vetCopia[i] = vetRecursos[i];
	}
	
	
    for(k=0; k<(l-1); k++){
        pivo = fabs(matrizCopia[k][k]); 
        posicao = k; 
        
        for(i=k+1; i<l; i++){
            if (fabs(matrizCopia[i][k]) > pivo){
				pivo = fabs(matrizCopia[i][k]); 
				posicao = i; 
            }
        }
		if (posicao != k){ 
		
			for (j=k; j<l; j++){ 
                aux = matrizCopia[k][j];
                matrizCopia[k][j] = matrizCopia[posicao][j];
                matrizCopia[posicao][j] = aux;
            }
            
            aux = vetCopia[k];
            vetCopia[k] = vetCopia[posicao];
            vetCopia[posicao] = aux;
        }
        
		//Metodo da eliminacao de Gauss (Triangularizacao)
        if (matrizCopia[k][k] != 0.0){ 
            for (i=k+1; i<l; i++){
				m = matrizCopia[i][k]/matrizCopia[k][k]; 
				matrizCopia[i][k] = 0.0;
				for(j=k+1; j<l; j++){ 
                	matrizCopia[i][j] = matrizCopia[i][j] - (m  * matrizCopia[k][j]);
                }
                vetCopia[i] = vetCopia[i] - m * vetCopia[k];
            }
        }
        else{
			erro = 1;
			printf("\n erro divisao por 0 \n");
			break;
		} 
	}
	
	if (erro == 1){ 
		liberarMemoriaMatriz(matrizCopia, l);
		free(vetCopia);
	}

	//Caso não haja erros na eliminação de Gauss, aplica o Metodo de resolucao de sistemas triangulares (Retrosubstituicao)
	
	if (matrizCopia[l-1][l-1] != 0.0){		
		vetSol[l-1] = vetCopia[l-1] / matrizCopia[l-1][l-1];	
	}
	    
	else{
		liberarMemoriaMatriz(matrizCopia,l);
		free(vetCopia); 
		printf("\n erro divisao por 0 \n");
	} 

	for (i=l-2; i>-1; i--){
		soma = 0.0;
		for (j=i+1; j<l; j++){
			soma = soma + matrizCopia[i][j]*vetSol[j];
		}

        if (matrizCopia[i][i] != 0){
			vetSol[i] = (vetCopia[i] - soma)/matrizCopia[i][i];
		}
        else{
			liberarMemoriaMatriz(matrizCopia, l);
			free(vetCopia);
			printf("\n erro divisao por 0 \n");
		} 
    }
	
	liberarMemoriaMatriz(matrizCopia,l);
	free(vetCopia);
}


//funcao que escolhe as matrizes Basica, Nao basica, custosBasico e custoNbasico da fase 2
void escolher_matrizBN(double **matrizA, double **matrizB, double **matrizN, int l, int c, double *vetCustos, double *vetCustoBasico, double *vetCustoNBasico, int escolha){
	
	int i, j, k;
	
	for(i = 0; i < l; i++){
		for(j = c-l, k = 0; j < c && k < c - (c-l); j++, k++){
			matrizB[i][k] = matrizA[i][j];
		}
	}
	
	printf("\nMatriz B:\n");
	

	for(i = 0; i < l; i++){
		for(j = 0; j < c - (c-l); j++){
			printf("%.2lf ", matrizB[i][j]);
		}
		printf("\n");
	}

	for(i = 0; i < l; i++){
		for(j = 0; j < c-l; j++){
			matrizN[i][j] = matrizA[i][j];
		}
	}
	
	printf("\nMatriz N:\n");
	
	for(i = 0; i < l; i++){
		for(j = 0; j < c-l; j++){
			printf("%.2lf ", matrizN[i][j]);
		}
		printf("\n");
	}

	for(i = c-l, k = 0; i < c && k < c - (c - l); i++, k++){
		vetCustoBasico[k] = vetCustos[i];
	}

	printf("\nVetor custo basico:\n");

	for(i = 0; i < l; i++){
		printf("%.2lf ", vetCustoBasico[i]);
	}

	for(i = 0; i < c-l ; i++){
		vetCustoNBasico[i] = vetCustos[i];
	}
	
	printf("\n\nVetor custo nao basico:\n");

	for(i = 0; i < c-l ; i++){
		printf("%.2lf ", vetCustoNBasico[i]);
	}

	printf("\n");

}


//funcao que percorre as linhas e colunas da matriz para determinar por qual fase comecamos o algoritmo
int percorre_matriz(double **matrizA, int l, int c){

	int linha, coluna, teste = 0, i;

	linha = 0;
	coluna = c - l;


	for(i = 0; i < l; i++){
		if(matrizA[linha][coluna] == 1){
			teste = 1;
			linha++;
			coluna++;
		}
		else{
			teste = 0;
			printf("Comecar algoritmo pela Fase 1:\n");
			return teste;
		}
	}

	return teste;
}


//impressao dos numeros das entradas inseridas formatado
void imprimir_numero_formatado(double num){

	if(num == 1){
		printf("+");
	}
	else if(num == -1){
		printf("-");
	}
	else if(num > 0){
		printf("+%.1f", num);
	}
	
	else if(num == 0){
		printf("+0");
	}

	else if(num < 0){
		printf("-%.1f", -num);
	}
	else{
		printf("+%.1f", num);
	}

}


//impressao de todas as entradas inseridas formatada
void printar_entrada_formatada(int linhas, int colunas, double *vetCustos, double *vetRecursos, double **matrizA){
	
	int i, j;

	printf("\nMin f(x) = ");

	for(i = 0; i < colunas; i++){
				
		imprimir_numero_formatado(vetCustos[i]);
	
		printf("X%d ", i+1);

	}

	printf("\n\n");
		
	printf("Sujeito A:\n");	

	for(i = 0; i < linhas; i++){
		for(j = 0; j < colunas; j++){

			imprimir_numero_formatado(matrizA[i][j]);
	
			printf("X%d ", j+1);
		}
		
		printf(" = %.1f\n", vetRecursos[i]);
	}

	for(i = 0; i < colunas; i++){
		if(i+1 < colunas){
			printf("X%d, ", i+1);
		}
		else{
			printf("X%d >= 0", i+1);
		}
	}
	
	
	
	printf("\n\n");
}


/*funcao main do programa*/


int main(){

	/*variaveis utilizadas no algoritmo*/ 
	int linhas, colunas, colunasR, identidade = 0, retorno = 0, testeNeg = 0, indice = 0, *trocas, contador = 0, controle = 0, escolha = 0, infactivel = 0;
	double *vetCustos, *vetCustosF1, *vetCustoBasico, *vetCustoNBasico, *vetRecursos, *vetSol, *vetLambda, *vetY, *vetColSel;
	double **matrizA, **matrizAF1, **matrizB, **matrizN;
	double menor = 0, resultado = 0;
	int i, j, testePos = 0, multiplas = 0, custozerado = 0;

	printf("\nObs.: As entradas devem ser digitadas na forma padrao, todas as variaveis sao positivas!\n\n");

	printf("Digite o numero de restricoes m e o numero de variaveis n separados por espaco:\n");

	scanf("%d %d", &linhas, &colunas);

	vetCustos = malloc(colunas * sizeof(double));
	vetRecursos = malloc(linhas * sizeof(double));
	vetSol = calloc(linhas, sizeof(double));
	vetLambda = calloc(linhas, sizeof(double));
	vetY = calloc(linhas, sizeof(double));
	vetColSel = calloc(linhas, sizeof(double));

	printf("Digite os valores c do vetor de custo:\n");	

	for(i = 0; i < colunas; i++){
		scanf("%lf", &vetCustos[i]);
	}
	
	printf("Digite os valores b do vetor de recursos:\n");	

	for(i = 0; i < linhas; i++){
		scanf("%lf", &vetRecursos[i]);
	}

	
	matrizA = malloc(linhas * sizeof(double));

	for(i = 0; i < linhas; i++){
		matrizA[i] = malloc(colunas * sizeof(double));
	}

	printf("Digite os valores da matriz A:\n");	
	
	for(i = 0; i < linhas; i++){
		for(j = 0; j < colunas; j++){
			scanf("%lf", &matrizA[i][j]);
		}
	}

	printf("\n\t**********Entrada formatada**********\t\n");

	printar_entrada_formatada(linhas, colunas, vetCustos, vetRecursos, matrizA);

	identidade = percorre_matriz(matrizA, linhas, colunas);

	if(identidade == 0){
		
		/*fase 1*/

		colunasR = colunas + linhas;

		vetCustosF1 = malloc(colunasR * sizeof(double));

		for(i = 0; i < colunasR; i++){
			if(i < colunas){
				vetCustosF1[i] = 0;
			}
			else{
				vetCustosF1[i] = 1;
			}
		}

		trocas = malloc(colunasR * sizeof(int));

		for(i = 0; i < colunasR; i++){
			trocas[i] = i;
		}

		matrizB = malloc(linhas * sizeof(double));
	
		for(i = 0; i < linhas; i++){
			matrizB[i] = malloc(linhas * sizeof(double));
		}
	
		matrizN = malloc(linhas * sizeof(double));
	
		for(i = 0; i < linhas; i++){
			matrizN[i] = malloc(colunas * sizeof(double));
		}

		
		vetCustoBasico = calloc((colunas-linhas), sizeof(double));
	
		vetCustoNBasico = malloc(colunas * sizeof(double));


		matrizAF1 = malloc(linhas*sizeof(double));
		for(i = 0; i < linhas; i++){
			matrizAF1[i] = malloc(colunasR * sizeof(double));
		}
		
		printf("\nMatriz A fase 1\n");
					
		for(i = 0; i < linhas; i++){
			for(j = 0; j < colunasR; j++){
				if(j < colunas){
					matrizAF1[i][j] = matrizA[i][j];
				}				
				else if(j >= colunas && j == (colunas + i)){
					matrizAF1[i][j] = 1;
				}
				else if(j >= colunas && j != (colunas + i)){
					matrizAF1[i][j] = 0;
				}

				printf("%.2lf ", matrizAF1[i][j]);
			}
			printf("\n");
		}

		printf("\nVetor custos fase 1\n");

		vetCustos = realloc(vetCustos, (colunasR) * sizeof(double));

		for(j = 0; j < colunasR; j++){		
			if(j >= colunas){
				vetCustos[j] = 0;
			}
			printf("%.2lf ", vetCustos[j]);
		}
		printf("\n");


		/*faz a fase 2 a partir do termino da fase 1 ou apenas a fase 2*/ 
		for(contador = 1; ; contador++){
			escolher_matrizBNfase1(matrizAF1, matrizB, matrizN, linhas, colunasR, vetCustosF1, vetCustoBasico, vetCustoNBasico);
			/*ira chamar todas as funcoes para fazer o metodo simplex
			primeiro chama a funcao void escolher_matrizBN*/

			resolver_sistema(matrizB, linhas, vetRecursos, vetSol); //metodo de gauss

			infactivel = 0;
			for(i = 0; i < linhas; i++){
				if(i == 0){
					infactivel = 0;
					printf("\nXB = %.2lf\n", vetSol[i]);
					if(vetSol[i] == 0){
						infactivel = 1;
					}			
				}
				else{
					printf("     %.2lf\n", vetSol[i]);
					if(vetSol[i] == 0){
						infactivel = 1;
					}				
				}
			}
			testePos = 0;	
			
			for(i = 0; i < linhas; i++){
				if(vetSol[i] == 0){
					testePos++;
				}
			}

			testeNeg = 0;
			
			for(i = 0; i < linhas; i++){
				if(vetSol[i] > 0){
					testeNeg++;
				}
			}
			
			if(testePos != 0){
				printf("\nSolucao basica degenerada\n\n");
				multiplas = 1;
			}
			else if(testeNeg == linhas){
				printf("\nSolucao basica nao degenerada\n\n");
				multiplas = 0;
			}

			resolver_sistema(transpor_matriz(matrizB, linhas), linhas, vetCustoBasico, vetLambda);
			
			for(i = 0; i < linhas; i++){
				if(i == 0){
					printf("\nVetor multiplicador simplex =\n%.2lf\n", vetLambda[i]);
				}
				else{
					printf("%.2lf\n", vetLambda[i]);
				}
			}
		
			retorno = custos_relativos(vetCustoNBasico, vetLambda, matrizN, linhas, colunasR, custozerado);
			
			if(retorno == -1){
				infactivel = 0;
				for(i = 0; i < linhas; i++){
					if(i == 0){
						printf("\nXB = %.2lf X%d\n", vetSol[i], trocas[colunasR-linhas+i]+1);
						if(vetSol[i] == 0){
							infactivel = 1;
						}
					}
					else{
						printf("     %.2lf X%d\n", vetSol[i], trocas[colunasR-linhas+i]+1);
						if(vetSol[i] == 0){
							infactivel = 1;
						}
					}
				}
				
				printf("\nf(x*) = %.2lf\n\n", multiplica(vetSol, vetCustoBasico, (colunasR - (colunasR - linhas))));
				if(multiplica(vetSol, vetCustoBasico, (colunasR - (colunasR - linhas))) == 0){
					escolha = 1;
					break;
				}
				else{
					printf("\nSolucao infactivel\n\n");
					return 0;
				}
			}
			else{
				
				for(i = 0; i < linhas; i++){
					vetColSel[i] = matrizN[i][retorno];
				}
				printf("\n");
				for(i = 0; i < linhas; i++){
					if(i == 0){
						printf("Coluna selecionada da matrizN %d = %.2lf\n", retorno+1, vetColSel[i]);
					}
					else{
						printf("                                  %.2lf\n", vetColSel[i]);
					}
				}
				
				resolver_sistema(matrizB, linhas, vetColSel, vetY);
				
				for(i = 0; i < linhas; i++){
					if(i == 0){
						printf("\nVet Y = %.2lf\n", vetY[i]); 
					}
					else{
						printf("        %.2lf\n", vetY[i]);
					}
				}
				
				testeNeg = 0;
				
				for(i = 0; i < linhas; i++){
					if(vetY[i] <= 0){
						testeNeg++;
					}
				}
				if(testeNeg == 0 && infactivel == 1){
					printf("\nSolucao infactivel\n\n");
					return 0;
				}
				controle = 0;
				for(i = 0; i < linhas; i++){
					if(vetY[i] > 0){
						resultado = vetSol[i]/vetY[i];
						printf("\nresultado XB%d/vetY%d %.2lf", i+1, i+1, resultado);
						if(controle == 0){
							menor = resultado;
							indice = i;
							controle = 1;
						}
						else if(resultado < menor){
							menor = resultado;
							indice = i;
						}
					}
				}
					
				printf("\n\nA coluna %d da matrizB sai da base!\n", indice+1);
					
				swapM(matrizAF1, retorno, (colunasR - linhas) + indice, linhas, trocas);
				swapV(vetCustosF1, retorno, (colunasR - linhas) + indice);
					
					
				for(i = 0; i < colunasR; i++){
					if(i == 0){
						printf("\nVet Custos trocado = %.2lf\n", vetCustosF1[i]);
					}
					else{
						printf("                     %.2lf\n", vetCustosF1[i]);
					}
				}
					
				for(i = 0; i < colunasR; i++){
					if(i == 0){
						printf("\ntroca das variaveis = %d\n", trocas[i]+1);
						
					}
					else{
						printf("                      %d\n", trocas[i]+1);
					}
				}
				
				printf("Matriz A trocada \n");
					
				for(i = 0; i < linhas; i++){
					for(j = 0; j < colunasR; j++){
						printf("%.2lf ", matrizAF1[i][j]);
					}
					printf("\n");
				}
				
				printf("\n\n**Fim da iteracao %d**", contador);
				printf("\n\nDigite qualquer numero inteiro para ir para a proxima iteracao:\n");

				scanf("%d", &i);
			}
		}
	}

	
	if(escolha == 0){
		printf("Comecar algoritmo pela Fase 2:\n");
	}
	else{
		printf("\n\n*************** Fase 2 ***************\n\n");
	}

	trocas = malloc(colunas * sizeof(int));

	for(i = 0; i < colunas; i++){
		trocas[i] = i;
	}

	matrizB = malloc(linhas * sizeof(double));
	
	for(i = 0; i < linhas; i++){
		matrizB[i] = malloc(linhas * sizeof(double));
	}

	matrizN = malloc (linhas * sizeof(double));

	for(i = 0; i < linhas; i++){
		matrizN[i] = malloc((colunas - linhas) * sizeof(double));
	}
	vetCustoBasico = calloc((colunas - linhas), sizeof(double));

	vetCustoNBasico = malloc((colunas - linhas) * sizeof(double));

	for(contador = 1; ;contador++){
		
		if(escolha == 1){
			escolher_matrizBNF2(matrizA, matrizB, matrizN, linhas, colunas, vetCustos, vetCustoBasico, vetCustoNBasico, escolha);
		}
		else{
			escolher_matrizBN(matrizA, matrizB, matrizN, linhas, colunas, vetCustos, vetCustoBasico, vetCustoNBasico, escolha);
		}
		
		/*ira chamar todas as funcoes p/ fazer o metodo simplex
		primeiro chama a funcao void escolher_matrizBN*/
		
		resolver_sistema(matrizB, linhas, vetRecursos, vetSol); //metodo de gauss
	
		for(i = 0; i < linhas; i++){			
			if(i == 0){
				infactivel = 0;
				printf("\nXB = %.2lf\n", vetSol[i]);
				if(vetSol[i] == 0){
					infactivel = 1;
				}			
			}
			else{
				printf("     %.2lf\n", vetSol[i]);
				if(vetSol[i] == 0){
					infactivel = 1;
				}			
			}
		}

		testePos = 0;	
			
		for(i = 0; i < linhas; i++){
			if(vetSol[i] == 0){
				testePos++;
			}
		}

		testeNeg = 0;
			
		for(i = 0; i < linhas; i++){
			if(vetSol[i] > 0){
				testeNeg++;
			}
		}

		if(testePos != 0){
			printf("\nSolucao basica degenerada\n\n");
		}
		else if(testeNeg == linhas){
			printf("\nSolucao basica nao degenerada\n\n");
		}

		resolver_sistema(transpor_matriz(matrizB, linhas), linhas, vetCustoBasico, vetLambda);
			
		for(i = 0; i < linhas; i++){
			if(i == 0){
				printf("\nVetor multiplicador simplex =\n%.2lf\n", vetLambda[i]);
			}
			else{
				printf("%.2lf\n", vetLambda[i]);
			}
		}
		
		retorno = custos_relativos(vetCustoNBasico, vetLambda, matrizN, linhas, colunas, custozerado);
			
		if(retorno == -1){
			
			if(escolha == 1){

				for(i = 0; i < linhas; i++){				
					if(vetSol[i] == 0){
						printf("\nSolucao otima degenerada\n\n");								
					}
				}
					printf("\nSolucao na iteracao atual e otima\n\n");

				for(i = 0; i < linhas; i++){
					if(i == 0){
						if(vetSol[i] == 0){
							printf("\nXB = 0.00 X%d\n", trocas[i]+1);
						}	
						else{									
							printf("\nXB = %.2lf X%d\n", vetSol[i], trocas[i]+1);
						}				
					}
					else{
						if(vetSol[i] == 0){
							printf("     0.00 X%d\n", trocas[i]+1);
						}	
						else{									
							printf("     %.2lf X%d\n", vetSol[i], trocas[i]+1);
						}	
					}
				}
			}

			else{
					for(i = 0; i < linhas; i++){				
						if(vetSol[i] == 0){
							printf("\nSolucao otima degenerada\n\n");								
						}
					}

					printf("\nSolucao na iteracao atual e otima\n\n");			
						
				for(i = 0; i < linhas; i++){
					if(i == 0){
						if(vetSol[i] == 0){
							printf("\nXB = 0.00 X%d\n", trocas[i+2]+1);
						}
						else{
							printf("\nXB = %.2lf X%d\n", vetSol[i], trocas[i+2]+1);
						}
					}
					else{
						if(vetSol[i] == 0){
							printf("     0.00 X%d\n", trocas[i+2]+1);
						}
						else{
							printf("     %.2lf X%d\n", vetSol[i], trocas[i+2]+1);
						}
					}
				}
			}
			
			printf("\nf(x*) = %.2lf\n\n", multiplica(vetSol, vetCustoBasico, (colunas - (colunas - linhas))));
			return 0;
		}

		else{
			
			for(i = 0; i < linhas; i++){
				vetColSel[i] = matrizN[i][retorno];
			}
			printf("\n");
			for(i = 0; i < linhas; i++){
				if(i == 0){
					printf("Coluna selecionada da matrizN %d = %.2lf\n", retorno+1, vetColSel[i]);
				}
				else{
					printf("                                  %.2lf\n", vetColSel[i]);
				}
			}
			
			resolver_sistema(matrizB, linhas, vetColSel, vetY);

			infactivel = 0;
			for(i = 0; i < linhas; i++){
				if(i == 0){
					printf("\nVet Y = %.2lf\n", vetY[i]);
					
				}
				else{
					printf("        %.2lf\n", vetY[i]);
				}
			}
			
			testeNeg = 0;
			
			for(i = 0; i < linhas; i++){
				if(vetY[i] <= 0){
					testeNeg++;
				}
			}


			
			if(testeNeg == linhas){
				printf("\nProblema nao tem solucao otima finita\n\n");
				return 0;
			}

			else if(testeNeg == 0 && infactivel == 1){
				printf("\nSolucao infactivel\n\n");
				return 0;
			}

			else{
				controle = 0;
				for(i = 0; i < linhas; i++){
					if(vetY[i] > 0){
						resultado = vetSol[i]/vetY[i];
						printf("\nresultado XB%d/vetY%d %.2lf", i+1, i+1, resultado);
						if(controle == 0){
							menor = resultado;
							indice = i;
							controle = 1;
						}
						else if(resultado < menor){
							menor = resultado;
							indice = i;
						}
					}
				}
			}	
				printf("\n\nA coluna %d da matrizB sai da base!\n", indice+1);
				
				if(escolha == 1){
					swapM(matrizA, colunas - (colunas - linhas) + retorno, indice, linhas, trocas);
					swapV(vetCustos, colunas - (colunas - linhas) + retorno, indice);
				}
				else{
					swapM(matrizA, retorno, (colunas - linhas) + indice, linhas, trocas);
					swapV(vetCustos, retorno, (colunas - linhas) + indice);
				}
				
				for(i = 0; i < colunas; i++){
					if(i == 0){
						printf("\nVet Custos trocado = %.2lf\n", vetCustos[i]);
					}
					else{
						printf("                     %.2lf\n", vetCustos[i]);
					}
				}
				
				for(i = 0; i < colunas; i++){
					if(i == 0){
						printf("\ntroca das variaveis = %d\n", trocas[i]+1);
						
					}
					else{
						printf("                      %d\n", trocas[i]+1);
					}
				}
			
				printf("Matriz A trocada \n");
				
				for(i = 0; i < linhas; i++){
					for(j = 0; j < colunas; j++){
						printf("%.2lf ", matrizA[i][j]);
					}
					printf("\n");
				}
					
				printf("\n\n**Fim da iteracao %d**", contador);
				printf("\n\nDigite qualquer numero inteiro para ir para a proxima iteracao:\n");

				scanf("%d", &i);
			}
	}
}
