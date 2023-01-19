import java.util.Arrays;
import java.util.HashMap;
import java.io.File;
import java.io.FileWriter;
import java.io.Writer;

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
    final static int POPULATION_SIZE = 1000;
    final static int MAX_GENERATIONS = 100;
    final static double CROSSOVER_RATE = 0.5;
    final static double MUTATION_RATE = 0.01;
    static boolean terminated; // Termination condition, flipped to true when termination condition is met
    static int generation; // Current generation

    // Population variables
    static int BITSTRING_LENGTH = 16; // Bitstring length
    static Population population = new Population(); // Population
    static Population children = new Population(); // Children
    static Individual champion; // Champion

    // fitness meausures
    static int maxFitness; // Max fitness
    static double avgFitness; // Average fitness
    static double identical; // Percentage of Identical individuals

    // File writer
    static Writer fileWriter;

    // main method
    public static void main(String[] args) {
        // Initialize variables
        generation = 0;
        terminated = false;

        try {
            fileWriter = new FileWriter(new File("output.txt"), false);
        } catch (Exception e) {
            System.out.println("Error: " + e);
        }

        // Print "one time header" <problem name> <population size> <bitstring genome length> <mutation rate> <crossover rate>
        writeToFile(PROBLEM_NAME + " " + POPULATION_SIZE + " " + BITSTRING_LENGTH + " " + MUTATION_RATE + " " + CROSSOVER_RATE + "\n");

        // Initialize population
        initializePopulation();

        // Print header for generation stats
        writeToFile("Generation\tChamp Fitness\tAvg Fitness\t\t%Identical");

        // Evaluate population
        evaluatePopulation();

        generation++;

        // While termination condition not met and max generations not reached
        while (!terminated) {
            // breed children by roulette wheel selection of parents and crossover
            breed();

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

        // Print champion after termination
        writeToFile("\nChampion genotype:\t" + champion.toString());
        writeToFile("Champion fitness:\t" + champion.getFitness());

        // close file writer
        try {
            fileWriter.close();
        } catch (Exception e) {
            System.out.println("Error: " + e);
        }
    }

    // Initialize population function
    public static void initializePopulation() {
        // Initialize population with uniform random values
        // build population with POPULATION_SIZE individuals,
        // each with a bit string of length BITSTRING_LENGTH, each bit randomly set to 0 or 1
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

    // calculate fitness of a single individual
    public static int calculateFitness(Individual individual) {
        // for maxones, fitness is the number of 1s in the bit string
        int fitness = 0;
        int[] bitString = individual.getBitString();
        for (int i = 0; i < BITSTRING_LENGTH; i++) {
            if (bitString[i] == 1) {
                fitness++;
            }
        }
        return fitness;
    }

    // Evaluate population function
    private static void evaluatePopulation() {
        // Evaluate population by calculating fitness of each individual and setting fitness variable
        Individual[] individuals = population.getIndividuals();
        for (int i = 0; i < POPULATION_SIZE; i++) {
            Individual individual = individuals[i];
            int fitness = calculateFitness(individual);
            individual.setFitness(fitness);
        }

        population.setIndividuals(individuals); // update population with fitness values

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
                champion = individual; // set champion to individual with max fitness
            }

            // average fitness total
            totalFitness += fitness;

            String bitStringString = individual.toString();
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
        writeToFile(generation + "\t\t\t" + maxFitness + "\t\t\t\t" + avgFitness + "\t\t\t\t" + identical);
    }

    // function to build roulette wheel (fitness proportionatal selection)
    private static Individual[] buildRouletteWheel(Individual[] individuals, int totalFitness) {
        Individual[] rouletteWheel = new Individual[totalFitness];
        
        int currentSlot = 0;

        // roulette wheen is a 0-totalFitness array of individuals
        // each individual is added to the roulette wheel a number of times equal to its fitness
        for (int i = 0; i < POPULATION_SIZE; i++) {
            Individual individual = individuals[i];
            int fitness = individual.getFitness();
            for (int j = 0; j < fitness; j++) {
                rouletteWheel[currentSlot] = individual;
                currentSlot++;
            }
        }

        return rouletteWheel;
    }

    // roulette select function
    private static Individual rouletteSelect(Individual[] rouletteWheel) {
        // select a random individual from the roulette wheel
        int randomIndex = (int) Math.floor(Math.random() * rouletteWheel.length);
        Individual selected = rouletteWheel[randomIndex];
        return selected;
    }

    // Select parents function
    private static void breed() {
        System.out.println("Selecting parents...");

        // get individuals from population
        Individual[] individuals = population.getIndividuals();

        // get total fitness
        int totalFitness = 0;
        for (int i = 0; i < POPULATION_SIZE; i++) {
            Individual individual = individuals[i];
            int fitness = individual.getFitness();
            totalFitness += fitness;
        }

        // build roulette wheel
        Individual[] rouletteWheel = buildRouletteWheel(individuals, totalFitness);

        Individual parent1 = null;
        Individual parent2 = null;

        // list of children
        Individual[] children = new Individual[POPULATION_SIZE];

        // Roulette wheel selection of parents
        for (int i = 0; i < POPULATION_SIZE; i++) {
            parent1 = rouletteSelect(rouletteWheel);
            parent2 = rouletteSelect(rouletteWheel);
            Individual child = crossover(parent1, parent2);
            children[i] = child;
        }

        // set children as population
        population.setIndividuals(children);
    }

    // Crossover function
    private static Individual crossover(Individual parent1, Individual parent2) {
        // TODO: Crossover (single point)

        return parent1;
    }

    // Mutate function
    private static void mutate() {
        // TODO: Mutate (bitwise with fixed mutation rate for each bit)
        System.out.println("Mutating...");
    }

    // Select survivors function
    private static void selectSurvivors() {
        // TODO: Select survivors (full replacement, all children replace all parents)
        System.out.println("Selecting survivors...");
    }

    // Check termination function
    private static void checkTermination() {
        // TODO: Check termination condition (max generations reached or population converged), print termination reason
        System.out.println("Checking termination...");
        if (generation >= MAX_GENERATIONS) {
            terminated = true;
        }
    }

    // write to file function
    private static void writeToFile(String line) {
        try {
            fileWriter.write(line + "\n");
        } catch (Exception e) {
            System.out.println("Error: " + e);
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

    // toString function
    public String toString() {
        return Arrays.toString(bitString);
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