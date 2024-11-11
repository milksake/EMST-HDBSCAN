# Fast Parallel Algorithms for Euclidean Minimum Spanning Tree (EMST) and Hierarchical Spatial Clustering

Este repositorio implementa los algoritmos paralelos rápidos para la construcción del Árbol de Expansión Mínima Euclidiana (EMST) y el Clustering Espacial Jerárquico, basados en el trabajo de Wang et al. (2021). El enfoque principal de este proyecto es la paralelización de algoritmos clásicos de geometría computacional para mejorar la eficiencia en el cálculo del EMST y el clustering jerárquico.

## Descripción del Proyecto

El algoritmo implementado en este repositorio es una adaptación de los métodos presentados en el artículo:

- **Título:** Fast Parallel Algorithms for Euclidean Minimum Spanning Tree and Hierarchical Spatial Clustering
- **Autores:** Yiqiu Wang, Shangdi Yu, Yan Gu, Julian Shun
- **Publicación:** Proceedings of the 2021 International Conference on Management of Data
- **DOI:** [10.1145/3448016.3457271](https://doi.org/10.1145/3448016.3457271)

Los algoritmos paralelizados incluyen:
1. **Cálculo del Árbol de Expansión Mínima Euclidiana (EMST)**
2. **Descomposición en Pares Bien Separados (WSPD)**
3. **Filtrado Geométrico para Kruskal**
4. **Clustering Espacial Jerárquico**

### Funciones Clave

- **Construcción Paralela del Árbol KD (KD-Tree):** Para organizar los puntos espaciales y facilitar el cálculo eficiente del EMST.
- **Cálculo de Distancias Euclidianas:** Determina la distancia entre nodos de manera paralela.
- **Kruskal Paralelizado:** Algoritmo para construir el EMST usando la técnica de Kruskal paralelizado.
- **WSPD (Well-Separated Pair Decomposition):** Algoritmo para encontrar pares de puntos bien separados en un conjunto de puntos.

## Instalación

### Requisitos

- **C++17** o superior
- **Compilador de C++** (como `g++` o `clang++`)
- **Hilos de ejecución (Multithreading)** soportado por tu plataforma.

### Pasos para la instalación

1. **Clonar el repositorio:**

   ```bash
   git clone https://github.com/tu-usuario/EMST-Parallel.git
   cd EMST-Parallel

2. **Compilar el proyecto:**

   Asegúrate de que tu compilador soporte C++17. Puedes usar `g++` para compilar el proyecto:

   ```bash
   g++ -o emst_parallel main.cpp

3. **Ejecutar el proyecto:**

    Una vez compilado, puedes ejecutar el proyecto de la siguiente manera:

   ```bash
   ./emst_parallel