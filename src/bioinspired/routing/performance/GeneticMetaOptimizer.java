package bioinspired.routing.performance;


import java.util.ArrayList;
import java.util.Random;

import bioinspired.routing.algorithm.RoutingParticleSwarm;
import bioinspired.routing.graph.GraphCreator;
import bioinspired.routing.routingTable.RoutingTable;

public class GeneticMetaOptimizer {

	/**
	 * SwarmSize: (Clerc2006_Page38 20-30) 10-50 iterations: 10-50 
	 * Initial Velocity: 0-20 
	 * W = (Clerc2006_Page41 0-1 smaller: premature convergence,larger: slow down convergence 0.7 or 0.8) 0-1 
	 * c1,c2 (Clerc2006_Page41 0-2, (0.7 1.47) (0.8 1.62)) 0-2 
	 * bestParticles = 10-SwarmSize
	 * 
	 * @param args
	 * 
	 */

	// constructor
	public GeneticMetaOptimizer(int populationSize, int generations,
			double mutationRate) {
		this.populationSize = populationSize;
		this.generations = generations;
		this.mutationRate = mutationRate;
		generateRandomInitialPopulation();
		geneticLoop();
		

	}

	public ArrayList<ArrayList<Double>> getPopulation() {
		return population;
	}
	
	public ArrayList<Double> getLowerFitnessPerGeneration() {
		return lowerFitnessPerGeneration;
	}
	
	public ArrayList<ArrayList<Double>> getLowerChromosomePerGeneration() {
		return lowerChromosomePerGeneration;
	}
	
	private void geneticLoop() {
		for (int i = 0; i < generations; i++) {
			evaluateChromosomes();
			selectParents();
			crossover();
			mutation();
			replacement();			
		}
		
	}

	// /********1**************
	private void generateRandomInitialPopulation() {
		for (int i = 0; i < populationSize; i++) {
			population.add(new ArrayList<Double>());
			double swarmSize = generateRandomIntegerNumber(10, 50);
			double iterations = generateRandomIntegerNumber(10, 50);
			double initialVelocity = generateRandomIntegerNumber(0, 20);
			double w = generateRandomDoubleNumber(1);
			double c1 = generateRandomDoubleNumber(2);
			double c2 = generateRandomDoubleNumber(2);
			double bestParticles = generateRandomIntegerNumber(10, 50);
			while (bestParticles > swarmSize) {
				bestParticles = generateRandomIntegerNumber(10, 50);
			}
			population.get(i).add(swarmSize);
			population.get(i).add(iterations);
			population.get(i).add(initialVelocity);
			population.get(i).add(w);
			population.get(i).add(c1);
			population.get(i).add(c2);
			population.get(i).add(bestParticles);
		}
	}

	private int generateRandomIntegerNumber(int min, int max) {
		Random rand = new Random();
		int num = rand.nextInt(max - min + 1) + min;
		return num;
	}

	private double generateRandomDoubleNumber(int max) {
		Random rand = new Random();
		double num = rand.nextDouble() * (max);
		return num;
	}

	// ***********2**************
	// fitness normalizado inversamente proporcional
	private void evaluateChromosomes() {
		double totalRawFitness = 0;
		for (int i = 0; i < populationSize; i++) {			
			double rawFitness = findFitness(population.get(i)); //original
			//double rawFitness = generateRandomDoubleNumber(50); //delete after test
			totalRawFitness += rawFitness;
			rawFitnessArray.add(rawFitness);
			System.out.print(i+" ");
			//System.out.println("Fitness "+rawFitness);			
		}
		findLowerFitness(rawFitnessArray);
		//System.out.println("Lower Chromosome "+lowerChromosomePerGeneration.toString());
		//System.out.println("Lower Fitness "+lowerFitnessPerGeneration.toString());
		double totalRawFitnessTemp = 0;
		for (int i = 0; i < rawFitnessArray.size(); i++) {
			double rawFitnessTemp = totalRawFitness / (rawFitnessArray.get(i));
			totalRawFitnessTemp += rawFitnessTemp;
			rawFitnessArrayTemp.add(rawFitnessTemp);
		}

		for (int i = 0; i < rawFitnessArrayTemp.size(); i++) {
			double fitness = (rawFitnessArrayTemp.get(i)) / totalRawFitnessTemp;
			FitnessArray.add(fitness);
		}
	}
	
	private void findLowerFitness(ArrayList<Double> fitnessArray) {
		double lowerFitness = fitnessArray.get(0);
		ArrayList<Double> lowerChromosome = population.get(0);
		for (int i = 0; i < fitnessArray.size(); i++) {
			double fitness = fitnessArray.get(i);
			if(fitness < lowerFitness) {
				lowerFitness = fitness;
				lowerChromosome = population.get(i);
			}			
		}		
		lowerFitnessPerGeneration.add(lowerFitness);
		lowerChromosomePerGeneration.add(lowerChromosome);		
	}
	
	private double findFitness(ArrayList<Double> chromosome) {
		
		int swarmSize = chromosome.get(0).intValue();
		int iterations = chromosome.get(1).intValue();
		int initialVelocity = chromosome.get(2).intValue();
		double w = chromosome.get(3);
		double c1 = chromosome.get(4);
		double c2 = chromosome.get(5);
		int bestParticles = chromosome.get(6).intValue();

		GraphCreator g = new GraphCreator("data/topology20NodesFull.txt", 
				"data/CostMatrix20Nodes(1-20).txt", "data/delayMatrix20Nodes(1-5).txt");		
		RoutingParticleSwarm rp = new RoutingParticleSwarm(swarmSize,
				iterations, w, c1, c2, initialVelocity, bestParticles, 15,
				5);
		RoutingTable rTable = new RoutingTable(g.getGraph());
		rTable.findKShortestPaths(0, 17, 200);
        rp.addRoutingTable(rTable.getRoutingTable());
        rTable.findKShortestPaths(0, 18, 200);
        rp.addRoutingTable(rTable.getRoutingTable());
        rTable.findKShortestPaths(0, 19, 200);
        rp.addRoutingTable(rTable.getRoutingTable());

		rp.addCostMatrix(g.getCostMatrix());
		rp.addDelayMatrix(g.getDelayMatrix());

		rp.generateInitialSwarm();
		rp.performParticleRouting();
		
		double rawFitness = rp.getMinFitness().get(rp.getMinFitness().size() - 1);
		
		return rawFitness;
		
	}

	// ********3*************
	// proportional selection scheme srinivas1994

	private void selectParents() {
		selectedParents = 0;
		while (true) {
			for (int i = 0; i < FitnessArray.size(); i++) {
				double FitnessInverseProportional = FitnessArray.get(i);
				Random rand = new Random();
				double randNumber = rand.nextDouble();
				// System.out.println("parent "+i);
				if (randNumber < FitnessInverseProportional) {
					parents.add(new ArrayList<Double>());
					parents.get(selectedParents).addAll(population.get(i));
					fitnessParents.add(FitnessArray.get(i));
					selectedParents++;
				}
				if (selectedParents == (populationSize / 2)) {
					break;
				}
				// System.out.println("selected "+selectedParents);
			}
			if (selectedParents == (populationSize / 2)) {
				break;
			}
		}
		//System.out.println("Parents " + " " + parents.size() + " "+ parents.toString());
	}

	// ***********4************
	// crossover one random point, random select a pair of parents from parents
	// set
	private void crossover() {
		for (int i = 0; i < populationSize / 2; i++) {
			int fatherID = generateRandomIntegerNumber(0, (parents.size() - 1));
			int motherID = generateRandomIntegerNumber(0, (parents.size() - 1));

			ArrayList<Double> father = parents.get(fatherID);
			ArrayList<Double> mother = parents.get(motherID);
			//System.out.println("Father " + father.toString());
			//System.out.println("mother " + mother.toString());

			int chromosomeLenght = father.size() - 1;
			int randomCrossoverPoint = generateRandomIntegerNumber(0,
					chromosomeLenght);
			//System.out.println(randomCrossoverPoint);
			// Crossover from father to mother
			sons1.add(new ArrayList<Double>());
			for (int j = 0; j <= randomCrossoverPoint; j++) {
				sons1.get(i).add(father.get(j));
			}
			for (int j = (randomCrossoverPoint + 1); j <= chromosomeLenght; j++) {
				sons1.get(i).add(mother.get(j));
			}

			// Crossover from mother to father
			sons2.add(new ArrayList<Double>());
			for (int j = 0; j <= randomCrossoverPoint; j++) {
				sons2.get(i).add(mother.get(j));
			}

			for (int j = (randomCrossoverPoint + 1); j <= chromosomeLenght; j++) {
				sons2.get(i).add(father.get(j));
			}
		}
		//System.out.println("sons1 " + sons1.toString());
		//System.out.println("sons2 " + sons2.toString());
		parentsAndSons.addAll(parents);
		parentsAndSons.addAll(sons1);
		parentsAndSons.addAll(sons2);

		//System.out.println("P and S " + " " + parentsAndSons.size() + "\t \t" + parentsAndSons.toString());
	}

	// ******5****

	// mutation: parents and sons pass through mutation

	private void mutation() {
		for (int i = 0; i < parentsAndSons.size(); i++) {
			for (int j = 0; j < parentsAndSons.get(0).size(); j++) {
				Random rand = new Random();
				double randNumber = rand.nextDouble();
				if (randNumber < mutationRate) {
					if (j == 0) {
						double swarmSize = generateRandomIntegerNumber(10, 50);
						parentsAndSons.get(i).set(j, swarmSize);

					} else if (j == 1) {
						double iterations = generateRandomIntegerNumber(10, 50);
						parentsAndSons.get(i).set(j, iterations);

					} else if (j == 2) {
						double initialVelocity = generateRandomIntegerNumber(0,
								20);
						parentsAndSons.get(i).set(j, initialVelocity);

					} else if (j == 3) {
						double w = generateRandomDoubleNumber(1);
						parentsAndSons.get(i).set(j, w);
					} else if (j == 4) {
						double c1 = generateRandomDoubleNumber(2);
						parentsAndSons.get(i).set(j, c1);
					} else if (j == 5) {
						double c2 = generateRandomDoubleNumber(2);
						parentsAndSons.get(i).set(j, c2);
					} else if (j == 6) {
						double bestParticles = generateRandomIntegerNumber(10,
								50);
						parentsAndSons.get(i).set(j, bestParticles);
					}
				}
			}
		}
		//System.out.println("P and S after mutation \t" + parents.toString());
	}
	
	private void replacement() {
		System.out.println("Replacement ");
		for (int i = 0; i < parentsAndSons.size(); i++) { // find the
																// fitness after
																// mutation
			fitnessAfterMutation.add(findFitness(parentsAndSons.get(i))); //original
			//fitnessAfterMutation.add(generateRandomDoubleNumber(50)); //delete after test
			
			System.out.print(i+" ");
		}
		for (int i = 0; i < populationSize; i++) { 
			double lowerFitness = fitnessAfterMutation.get(0);
			ArrayList<Double> lowerChromosome = parentsAndSons.get(0);
			int position = 0;

			for (int j = 0; j < fitnessAfterMutation.size(); j++) {
				double fitness = fitnessAfterMutation.get(j);
				if (fitness < lowerFitness) {
					lowerFitness = fitness;	
					lowerChromosome = parentsAndSons.get(j);
					position = j;
				}
			}
			populationAfterMutation.add(lowerChromosome);
			fitnessAfterMutation.set(position, Double.POSITIVE_INFINITY);
		}
		//System.out.println("population after mutation"+" "+populationAfterMutation.size()+" "+populationAfterMutation.toString());
		population.clear();
		population.addAll(populationAfterMutation);
		rawFitnessArray.clear();
		rawFitnessArrayTemp.clear();
		FitnessArray.clear();
		parents.clear();
		fitnessParents.clear();
		sons1.clear();
		sons2.clear();
		parentsAndSons.clear();
		fitnessAfterMutation.clear();
		populationAfterMutation.clear();
	}

	// instance variables
	private int populationSize;
	private int generations;
	private double mutationRate;
	private int selectedParents = 0;
	private ArrayList<ArrayList<Double>> population = new ArrayList<ArrayList<Double>>();
	private ArrayList<Double> rawFitnessArray = new ArrayList<Double>();
	private ArrayList<Double> rawFitnessArrayTemp = new ArrayList<Double>();
	private ArrayList<Double> FitnessArray = new ArrayList<Double>();
	private ArrayList<Double> fitnessParents = new ArrayList<Double>();
	private ArrayList<ArrayList<Double>> parents = new ArrayList<ArrayList<Double>>();
	private ArrayList<ArrayList<Double>> sons1 = new ArrayList<ArrayList<Double>>();
	private ArrayList<ArrayList<Double>> sons2 = new ArrayList<ArrayList<Double>>();
	private ArrayList<ArrayList<Double>> parentsAndSons = new ArrayList<ArrayList<Double>>();
	private ArrayList<ArrayList<Double>> populationAfterMutation = new ArrayList<ArrayList<Double>>();
	private ArrayList<Double> fitnessAfterMutation = new ArrayList<Double>();
	private ArrayList<Double> lowerFitnessPerGeneration = new ArrayList<Double>();
	private ArrayList<ArrayList<Double>> lowerChromosomePerGeneration = new ArrayList<ArrayList<Double>>(); 

}
