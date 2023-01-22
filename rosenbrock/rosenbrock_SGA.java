package rosenbrock;

/*
 * Jackson Hacker
 * CS 5146: Evolutionary Computing
 * Simple Genetic Algorithm
 * 1/17/2023
 */


// Import statements
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.io.File;
import java.io.FileWriter;
import java.io.Writer;
import java.util.ArrayList;

// Class definition
class rosenbrock_SGA {
    // final Static variables
    final static String PROBLEM_NAME = "rosenbrock";
    final static int POPULATION_SIZE = 500;
    final static int MAX_GENERATIONS = 2000;
    final static double CROSSOVER_RATE = 0.7;
    final static double MUTATION_RATE = 0.05;
    final static int GENOME_LENGTH = 16; // Bitstring length
    final static int CONVERGENCE_THRESHOLD = 5; // Convergence threshold (number of generations with similar average fitness to terminate)
    final static double X_BOUND = 5; // x values can be between -N to N
    final static double Y_BOUND = 5; // y values can be between -N to N

    static boolean terminated; // Termination condition, flipped to true when termination condition is met
    static int generation; // Current generation

    // Population variables
    static Population population = new Population(); // Population
    static Population children = new Population(); // Children
    static Individual champion; // Champion

    // fitness meausures
    static double maxFitness; // Max fitness
    static double avgFitness; // Average fitness
    static double identical; // Percentage of Identical individuals

    static double[] avgHistory = new double[CONVERGENCE_THRESHOLD]; // History of average fitness over past N generations
    
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
        writeToFile(PROBLEM_NAME + " " + POPULATION_SIZE + " " + GENOME_LENGTH + " " + MUTATION_RATE + " " + CROSSOVER_RATE + "\n");

        // Initialize population
        initializePopulation();

        // Print header for generation stats
        writeToFile("Generation\tChamp Fitness\t\t\t\t\tAvg Fitness\t\t\t%Identical");

        // Evaluate population
        evaluatePopulation();

        generation++;

        // While termination condition not met and max generations not reached
        while (!terminated) {
            // breed children by roulette wheel selection of parents and crossover
            breedAndReplace(); // also replaces population with children

            // Mutate
            mutate();

            // Evaluate population
            evaluatePopulation();


            // Check termination condition (flip terminated to true if met)
            checkTermination();

            // Increment iterations
            generation++;
        }

        // Print champion after termination
        writeToFile("\nChampion genotype:\t" + champion.toString());
        writeToFile("Champion fitness:\t" + champion.getFitness());
        writeToFile("Champion f(x, y):\t" + -1 * champion.getFitness());
        writeToFile("Champion x value:\t" + getX(champion.getBitString()));
        writeToFile("Champion y value:\t" + getY(champion.getBitString()));

        // close file writer
        try {
            fileWriter.close();
        } catch (Exception e) {
            System.out.println("Error: " + e);
        }

        // Print location of output file
        System.out.println("Output file written to: " + new File("output.txt").getAbsolutePath());
    }

    // Initialize population function
    public static void initializePopulation() {
        // Initialize population with uniform random values
        // build population with POPULATION_SIZE individuals,
        // each with a bit string of length BITSTRING_LENGTH, each bit randomly set to 0 or 1
        Individual[] individuals = new Individual[POPULATION_SIZE];
        for (int i = 0; i < POPULATION_SIZE; i++) {
            Individual individual = new Individual();
            int[] bitString = new int[GENOME_LENGTH];
            for (int j = 0; j < GENOME_LENGTH; j++) {
                bitString[j] = (int) Math.round(Math.random());
            }
            individual.setBitString(bitString);
            individuals[i] = individual;
        }
        population.setIndividuals(individuals);
    }

    // binary to decimal function
    public static int binaryToDecimal(int[] binary) {
        int decimal = 0;
        for (int i = 0; i < binary.length; i++) {
            decimal += binary[i] * Math.pow(2, binary.length - i - 1);
        }
        return decimal;
    }

    // rosnebrock function
    public static double rosenbrock(double x, double y) {
        // f(x,y) = (a-x)^2 + b(y-x^2)^2
        // a = 1, b = 100

        double a = 1;
        double b = 100;

        return Math.pow(a - x, 2) + b * Math.pow(y - Math.pow(x, 2), 2);
    }

    // get x value from bitstring
    public static double getX(int[] bitString) {
        // get x value from bitstring
        int[] x = Arrays.copyOfRange(bitString, 0, GENOME_LENGTH / 2);
        int xDecimal = binaryToDecimal(x);
        // scale x value to be between -X_BOUND and X_BOUND
        return (-1 * X_BOUND) + (xDecimal * 10.0 / (Math.pow(2, GENOME_LENGTH / 2) - 1));
    }

    // get y value from bitstring
    public static double getY(int[] bitString) {
        // get y value from bitstring
        int[] y = Arrays.copyOfRange(bitString, GENOME_LENGTH / 2, GENOME_LENGTH);
        int yDecimal = binaryToDecimal(y);
        return (-1 * Y_BOUND) + (yDecimal * 10.0 / (Math.pow(2, GENOME_LENGTH / 2) - 1));
    }

    // calculate fitness of a single individual
    public static double calculateFitness(Individual individual) {
        // for minimizing rosenbrock, the fitness is the negative of the rosenbrock function (maximize the negative to minimize)
        // from our representation, the first half of the bitstring is the x value, and the second half is the y value
        // the x and y values are converted to a decimal value between -X_BOUND and X_BOUND, -Y_BOUND and Y_BOUND
        // the rosenbrock function is then calculated with the x and y values

        // get bitstring
        int[] bitString = individual.getBitString();

        double xValue = getX(bitString);
        double yValue = getY(bitString);

        // calculate rosenbrock function
        double fitness = -1*rosenbrock(xValue, yValue);

        return fitness;
    }

    // Evaluate population function
    private static void evaluatePopulation() {
        // Evaluate population by calculating fitness of each individual and setting fitness variable
        Individual[] individuals = population.getIndividuals();
        for (int i = 0; i < POPULATION_SIZE; i++) {
            Individual individual = individuals[i];
            double fitness = calculateFitness(individual);
            individual.setFitness(fitness);
        }

        population.setIndividuals(individuals); // update population with fitness values

        // total fitness
        int totalFitness = 0;
        // indistinct individuals
        int indistinct = 0;

        // Hash map of distinct individuals
        HashMap<String, Integer> distincts = new HashMap<String, Integer>();

        // set max fitness to first individual's fitness
        maxFitness = individuals[0].getFitness();

        // Calculate max fitness, average fitness, and percentage of identical individuals
        for (int i = 0; i < POPULATION_SIZE; i++) {
            Individual individual = individuals[i];
            double fitness = individual.getFitness();

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

        // update average history
        for (int i = 0; i < avgHistory.length - 1; i++) {
            avgHistory[i] = avgHistory[i + 1];
        }
        avgHistory[avgHistory.length - 1] = avgFitness;

        // Print max fitness, average fitness, and percentage of identical individuals
        // <gen number> <fitness score of champion> <average fitness score> <percent of genomes in pop that are identical>
        writeToFile(generation + "\t\t\t" + maxFitness + "\t\t\t\t" + avgFitness + "\t\t\t\t" + identical);
    }

     // function to build roulette wheel (fitness proportionatal selection)
    private static List<Individual> buildRouletteWheel(Individual[] individuals, double totalFitness) {
        // for each individual, add it to the roulette wheel a proportional amount of times according to its fitness
        // the more fit an individual is, the more times it will be added to the roulette wheel
        List <Individual> rouletteWheel = new ArrayList<Individual>();
        // sort the list (ascending)
        Arrays.sort(individuals, new Comparator<Individual>() {
            @Override
            public int compare(Individual i1, Individual i2) {
                return Double.compare(i1.getFitness(), i2.getFitness());
            }
        });

        // add individuals to roulette wheel, proportional to fitness (smallest fitness added zero times, largest fitness added POPULATION_SIZE-1 times)
        for (int i = 0; i < POPULATION_SIZE; i++) {
            Individual individual = individuals[i];
            int numTimes = i;
            for (int j = 0; j < numTimes; j++) {
                rouletteWheel.add(individual);
            }
        }

        return rouletteWheel;
    }

    // roulette select function
    private static Individual rouletteSelect(List<Individual> rouletteWheel) {
        // select a random individual from the roulette wheel
        int randomIndex = (int) Math.floor(Math.random() * rouletteWheel.size());
        Individual selected = rouletteWheel.get(randomIndex);
        return selected;
    }

    // returns children from roulette wheel selection of parents and crossover
    // also replaces population with children (full replacement)
    private static void breedAndReplace() {
        // get individuals from population
        Individual[] individuals = population.getIndividuals();

        // get total fitness
        double totalFitness = 0;
        for (int i = 0; i < POPULATION_SIZE; i++) {
            Individual individual = individuals[i];
            double fitness = individual.getFitness();
            totalFitness += fitness;
        }

        // build roulette wheel
        List<Individual> rouletteWheel = buildRouletteWheel(individuals, totalFitness);

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
        // generate random number between 0 and 1, if less than crossover rate, crossover
        // if not, return parent1 (clone)
        double random = Math.random();
        if (random < CROSSOVER_RATE) {
            int[] parent1BitString = parent1.getBitString();
            int[] parent2BitString = parent2.getBitString();

            // generate random crossover point
            int crossoverPoint = (int) Math.floor(Math.random() * GENOME_LENGTH);

            // create child bit string
            int[] childBitString = new int[GENOME_LENGTH];

            // copy first part of parent1 bit string to child bit string
            for (int i = 0; i < crossoverPoint; i++) {
                childBitString[i] = parent1BitString[i];
            }

            // copy second part of parent2 bit string to child bit string
            for (int i = crossoverPoint; i < GENOME_LENGTH; i++) {
                childBitString[i] = parent2BitString[i];
            }

            // create child individual
            Individual child = new Individual();
            child.setBitString(childBitString);
            return child;
        }

        return parent1;
    }

    // Mutate function
    private static void mutate() {
        // mutate all individuals in population, flipping each bit with a probability of MUTATION_RATE
        Individual[] individuals = population.getIndividuals();
        for (int i = 0; i < POPULATION_SIZE; i++) {
            Individual individual = individuals[i];
            int[] bitString = individual.getBitString();
            for (int j = 0; j < GENOME_LENGTH; j++) {
                double random = Math.random();
                if (random < MUTATION_RATE) { // flip if random number is less than mutation rate
                    if (bitString[j] == 0) {
                        bitString[j] = 1;
                    } else {
                        bitString[j] = 0;
                    }
                }
            }
            individual.setBitString(bitString);
        }
    }

    // Check termination function
    private static void checkTermination() {
        // check if max generations reached
        if (generation >= MAX_GENERATIONS) {
            terminated = true;

            // write to file termination reason
            writeToFile("TERMINATED AT GENERATION " + generation + "\n");
            writeToFile("Termination reason: Max generations reached\n\n");
        }

        // check if population converged 
        //(average fitness of population has not changed by more than 1% in the last CONVERGENCE_THRESHOLD generations)
        if (generation >= CONVERGENCE_THRESHOLD) {
            // if percent change in average fitness is less than 1% between each pair of generations in the last CONVERGENCE_THRESHOLD generations, terminate
            double percentChange = 0;
            boolean converged = true;
            for (int i = 0; i < CONVERGENCE_THRESHOLD - 1; i++) {
                double fitness1 = avgHistory[i];
                double fitness2 = avgHistory[i + 1];
                percentChange = Math.abs((fitness2 - fitness1) / fitness1);
                if (percentChange > 0.01) {
                    converged = false;
                    break;
                }
            }
            if (converged) {
                terminated = true;

                // write to file termination reason
                writeToFile("TERMINATED AT GENERATION " + generation + "\n");
                writeToFile("Termination reason: Average fitness converged\n\n");
            }
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
    double fitness; // Fitness

    // Constructor
    public Individual() {
        // Initialize bit string
        bitString = new int[rosenbrock_SGA.GENOME_LENGTH];

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
    public double getFitness() {
        return fitness;
    }

    // Set fitness function
    public void setFitness(double fitness) {
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
        individuals = new Individual[rosenbrock_SGA.POPULATION_SIZE];
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