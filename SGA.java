import java.util.Arrays;
import java.util.HashMap;

/*
 * Jackson Hacker
 * CS 5146: Evolutionary Computing
 * Simple Genetic Algorithm
 * 1/17/2023
 */


// Import statements


// Class definition
class SGA {
    // Static variables
    final static String PROBLEM_NAME = "maxones";
    final static int POPULATION_SIZE = 100;
    final static int MAX_GENERATIONS = 1;
    final static double CROSSOVER_RATE = 0.5;
    final static double MUTATION_RATE = 0.01;
    static boolean terminated; // Termination condition, flipped to true when termination condition is met
    static int generation; // Current generation

    // Population variables
    static int BITSTRING_LENGTH = 8; // Bitstring length
    static Population population; // Population
    static Population parents; // Parents
    static Population children; // Children

    // fitness meausures
    static int maxFitness; // Max fitness
    static double avgFitness; // Average fitness
    static double identical; // Percentage of Identical individuals

    // main method
    public static void main(String[] args) {
        // Initialize variables
        generation = 0;
        terminated = false;

        // Print "one time header" <problem name> <population size> <bitstring genome length> <mutation rate> <crossover rate>
        System.out.println(PROBLEM_NAME + " " + POPULATION_SIZE + " " + BITSTRING_LENGTH + " " + MUTATION_RATE + " " + CROSSOVER_RATE);

        // Initialize population
        initializePopulation();

        // Print header for generation stats
        System.out.println("Generation\tChamp Fitness\tAvg Fitness\t%Identical");

        // Evaluate population
        evaluatePopulation();

        generation++;

        // While termination condition not met and max generations not reached
        while (!terminated) {
            // Print current generation
            System.out.println("Generation " + generation);

            // Select parents
            selectParents();

            // Crossover
            crossover();

            // Mutate
            mutate();

            // Evaluate population
            evaluatePopulation();

            // Select survivors
            selectSurvivors();

            // Check termination condition (flip terminated to true if met)
            checkTermination();

            // Increment iterations
            generation++;
        }
    }

    // Initialize population function
    public static void initializePopulation() {
        // Initialize population with uniform random values
        // build population with POPULATION_SIZE individuals,
        // each with a bit string of length BITSTRING_LENGTH, each bit randomly set to 0 or 1
        population = new Population();
        Individual[] individuals = new Individual[POPULATION_SIZE];
        for (int i = 0; i < POPULATION_SIZE; i++) {
            Individual individual = new Individual();
            int[] bitString = new int[BITSTRING_LENGTH];
            for (int j = 0; j < BITSTRING_LENGTH; j++) {
                bitString[j] = (int) Math.round(Math.random());
            }
            individual.setBitString(bitString);
            individuals[i] = individual;
        }
        population.setIndividuals(individuals);
    }

    // Evaluate population function
    public static void evaluatePopulation() {
        // Evaluate population by calculating fitness of each individual and setting fitness variable
        // For max one problem, fitness is the number of 1s in the bit string
        Individual[] individuals = population.getIndividuals();
        for (int i = 0; i < POPULATION_SIZE; i++) {
            Individual individual = individuals[i];
            int[] bitString = individual.getBitString();
            int fitness = 0;
            for (int j = 0; j < BITSTRING_LENGTH; j++) {
                if (bitString[j] == 1) {
                    fitness++;
                }
            }
            individual.setFitness(fitness);
        }

        population.setIndividuals(individuals); // update population with new fitness values

        // total fitness
        int totalFitness = 0;
        // indistinct individuals
        int indistinct = 0;

        // Hash map of distinct individuals
        HashMap<String, Integer> distincts = new HashMap<String, Integer>();

        // Calculate max fitness, average fitness, and percentage of identical individuals
        for (int i = 0; i < POPULATION_SIZE; i++) {
            Individual individual = individuals[i];
            int fitness = individual.getFitness();

            // max fitness
            if (fitness > maxFitness) {
                maxFitness = fitness;
            }

            // average fitness total
            totalFitness += fitness;


            int[] bitString = individual.getBitString();
            String bitStringString = Arrays.toString(bitString);
            if (distincts.containsKey(bitStringString)) {
                // if individual is already in the hash map, increment the count (found duplicate)
                indistinct++;
            } else {
                // if individual is not in the hash map, add it to the hash map
                distincts.put(bitStringString, 1);
            } 
        }

        // average fitness final calculation
        avgFitness = (double) totalFitness / POPULATION_SIZE;

        // percentage of identical individuals
        identical = ((double) indistinct / POPULATION_SIZE) * 100;

        // Print max fitness, average fitness, and percentage of identical individuals
        // <gen number> <fitness score of champion> <average fitness score> <percent of genomes in pop that are identical>
        System.out.println(generation + "\t\t" + maxFitness + "\t\t" + avgFitness + "\t\t" + identical);
    }

    // Select parents function
    public static void selectParents() {
        // TODO: Select parents via roulette wheel selection (fitness proportionatal)
        System.out.println("Selecting parents...");
    }

    // Crossover function
    public static void crossover() {
        // TODO: Crossover (single point)
        System.out.println("Crossover...");
    }

    // Mutate function
    public static void mutate() {
        // TODO: Mutate (bitwise with fixed mutation rate for each bit)
        System.out.println("Mutating...");
    }

    // Select survivors function
    public static void selectSurvivors() {
        // TODO: Select survivors (full replacement, all children replace all parents)
        System.out.println("Selecting survivors...");
    }

    // Check termination function
    public static void checkTermination() {
        // TODO: Check termination condition (max generations reached)
        System.out.println("Checking termination...");
        if (generation >= MAX_GENERATIONS) {
            terminated = true;
        }
    }

}

// individual class (fixed length bit string)
class Individual {
    // Instance variables
    int[] bitString; // Bit string
    int fitness; // Fitness

    // Constructor
    public Individual() {
        // Initialize bit string
        bitString = new int[SGA.BITSTRING_LENGTH];

        // Initialize fitness
        fitness = 0;
    }

    // Get bit string function
    public int[] getBitString() {
        return bitString;
    }

    // Set bit string function
    public void setBitString(int[] bitString) {
        this.bitString = bitString;
    }

    // Get fitness function
    public int getFitness() {
        return fitness;
    }

    // Set fitness function
    public void setFitness(int fitness) {
        this.fitness = fitness;
    }
}

// population class (collection of individuals)
class Population {
    // Instance variables
    Individual[] individuals; // Individuals

    // Constructor
    public Population() {
        // Initialize individuals
        individuals = new Individual[SGA.POPULATION_SIZE];
    }

    // Get individuals function
    public Individual[] getIndividuals() {
        return individuals;
    }

    // Set individuals function
    public void setIndividuals(Individual[] individuals) {
        this.individuals = individuals;
    }
}