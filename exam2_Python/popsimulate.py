#_________________ Función 1
import numpy as np           #Biblioteca Numpy as np
def build_population(N, p): #Se define la función build_population(N, p) con N como tamaño y p probabilidad de obtener un alelo en concreto
    population = []        #Lista vacía
    for i in range(N):     #Bucle for para iterar N veces
        allele1 = "A"
        if np.random.rand() > p:
            allele1 = "a"
        allele2 = "A"
        if np.random.rand() > p:
            allele2 = "a"
        population.append((allele1, allele2))
    return population


#_________________ Función 2
def compute_frequencies(population):   #Se define la función con population para la población genética
    AA = population.count(("A", "A"))  
    Aa = population.count(("A", "a"))
    aA = population.count(("a", "A"))
    aa = population.count(("a", "a"))
    return({"AA": AA, "aa": aa, "Aa": Aa, "aA": aA})  #Se devuelve un diccionario con la fecuencia de los genotipos como principales ("AA", "aa", "Aa", "aA")

#_________________ Función 3
def reproduce_population(population): #Defino la función
    new_generation = []               #Lista vacía
    N = len(population)               #Tamaño de la población
    for i in range(N):                #Bucle for que itera N veces para generar N descendientes en nueva generación.
        dad = np.random.randint(N)    #Se selecciona aleatoriamente un padre y una madre seleccionados para ser los padres
        mom = np.random.randint(N)
        chr_mom = np.random.randint(2) #Selecciona aleatoriamente uno de los dos alelos de la madre. chr_mom tomará el valor de 0 o 1
        offspring = (population[mom][chr_mom], population[dad][1 - chr_mom]) #Combinación de Alelos para nuevo individuo
        new_generation.append(offspring)  #El nuevo individuo se genera a la lista new_generation
    return new_generation

#_________________ Función 4
def simulate_drift(N, p):   #Se define la función
    my_pop = build_population(N, p)   #Se llama a la función para crear una población inicial representada por my_pop.
    fixation = False                 #Fixation como fijador del alelo y num_generation como contador de generaciones
    num_generations = 0
    while fixation == False:         #Se inicia un bucle while donde se llama la funcion en genotype_counts
        genotype_counts = compute_frequencies(my_pop)#Se verifica si uno de los genotipos (AA o aa) ha alcanzado la fijación en la población
        if (genotype_counts["AA"] == N or genotype_counts["aa"] == N):   #Si se cumple esta condición se imprime un mensaje
            print("An allele reached fixation at generation", num_generations)  
            print("The genotype counts are")
            print(genotype_counts)     #Se muestra el recuento de genotipos
            fixation == True           #Se establece fixation en True para salir del bucle
            break
        my_pop = reproduce_population(my_pop)  #Se llama la funcion en my_pop
        num_generations += 1                   #Se incrementa el contador num_generations para seguir la cuenta de las generaciones
    return num_generations, genotype_counts  #Se devuelve el número total de generaciones y el recuento de genotipos en forma de tupla