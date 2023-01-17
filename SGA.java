/*
 * Jackson Hacker
 * CS 5146: Evolutionary Computing
 * Simple Genetic Algorithm
 * 1/17/2023
 */


// Import statements


// Class definition
class SGA {
    // Global variables
    final static int POPULATION_SIZE = 100;
    final static int MAX_GENERATIONS = 10;
    static boolean terminated; // Termination condition, flipped to true when termination condition is met
    static int generation; // Current generation

    // main method
    public static void main(String[] args) {
        // Initialize variables
        generation = 0;
        terminated = false;

        // Initialize population
        initializePopulation();

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
        // TODO: Initialize population with uniform random values
        System.out.println("Initializing population...");
    }

    // Evaluate population function
    public static void evaluatePopulation() {
        // TODO: Evaluate population
        System.out.println("Evaluating population...");
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