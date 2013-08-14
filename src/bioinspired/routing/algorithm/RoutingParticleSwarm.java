package bioinspired.routing.algorithm;

import java.util.ArrayList;

/**
 * This class is the core of the Multicast Routing Algorithm with QoS constrains
 * based on Particle Swarm Optimization. It uses the routing table built with
 * the K-Shortest Path algorithm in order create the particles and perform the
 * PSO algorithm to find the multicast tree which satisfies the QoS constrains
 * 
 * @author YEISON JULIAN CAMARGO BARAJAS
 */
public class RoutingParticleSwarm {

	/*
	 * Constructor
	 */
	public RoutingParticleSwarm(int swarmSize, int iterations) {
		i_swarmSize = swarmSize;
		i_iterations = iterations;
	}

	public RoutingParticleSwarm(int swarmSize, int iterations, double w,
			double c1, double c2, int initialVelocity, int bestParticles,
			double maxDelayOfMulticastTree, double maxDelayVariation) {
		i_swarmSize = swarmSize;
		i_iterations = iterations;
		this.w = w;
		this.c1 = c1;
		this.c2 = c2;
		this.initialVelocity = initialVelocity;
		this.bestParticles = bestParticles;
		delayMulticast = maxDelayOfMulticastTree;
		delayVariation = maxDelayVariation;
	}

	/*
	 * Getters and Setters
	 */
	/**
	 * The method returns the Nb number of best particles. This Nb set of
	 * particles is used to update the position of the particles based on the
	 * sigmod function using the elite approach.
	 * 
	 * @return An integer that represents the number of Nb set in the
	 *         constructor. If not value is passed to the constructor the
	 *         default value is 20.
	 */

	public int getEliteGroup() {
		return bestParticles;
	}

	/**
	 * This method allows to set the value Nb (number of best particles). The Nb
	 * set of particles is used to update the position of the particles based on
	 * the sigmod function using the elite approach.
	 * 
	 * @param bestParticles
	 *            An integer to tell the algorithm the quantity of the particles
	 *            that will be taken into account to be selected as the elite
	 *            group.
	 */

	public void setEliteGroup(int bestParticles) {
		this.bestParticles = bestParticles;
	}

	/**
	 * This method returns the cognitive constant set in the constructor. If not
	 * parameter is passed the default value is 1.49445. The cognitive constant
	 * represent the self-trust of the position found by the particle.
	 * 
	 * @return An integer value representing the cognitive constant
	 */

	public double getC1() {
		return c1;
	}

	/**
	 * This method allows to set the value of the cognitive constant. The
	 * cognitive constant represent the self-trust of the position found by the
	 * particle.
	 * 
	 * @param c1
	 *            An integer to set the cognitive constant. It must be greater
	 *            than zero. Normally c1 = c2
	 */

	public void setC1(double c1) {
		this.c1 = c1;
	}

	/**
	 * This method returns the social constant set in the constructor. If not
	 * parameter is passed the default value is 1.49445. The social constant
	 * represent the trust in the swarm of the position found by the swarm.
	 * 
	 * @return An integer value representing the social constant
	 */

	public double getC2() {
		return c2;
	}

	/**
	 * This method allows to set the value of the social constant. The social
	 * constant represent the trust in the swarm of the position found by the
	 * swarm.
	 * 
	 * @param c2
	 *            An integer to set the social constant. It must be greater than
	 *            zero. Normally c1 = c2
	 */

	public void setC2(double c2) {
		this.c2 = c2;
	}

	/**
	 * This method allows to get the number of iterations set in the
	 * constructor. The iterations represent the times the algorithm will be
	 * executed.
	 * 
	 * @return An integer that represents the number of iterations
	 */

	public int getNumberOfIterations() {
		return i_iterations;
	}

	/**
	 * This method allows to set the number of iterations. The iterations
	 * represent the times the algorithm will be executed.
	 * 
	 * @param iterations
	 *            An integer that represents the number of iterations
	 */

	public void setNumberOfIterations(int iterations) {
		i_iterations = iterations;
	}

	/**
	 * The method returns the size of the swarm. The size of the swarm dictates
	 * the number of particles that will be included in the swarm.
	 * 
	 * @return The size of the swarm.
	 */

	public int getSwarmSize() {
		return i_swarmSize;
	}

	/**
	 * The method allows to set the size of the swarm.
	 * 
	 * @param swarmSize
	 *            The size of the swarm.
	 */

	public void setSwarmSize(int swarmSize) {
		i_swarmSize = swarmSize;
	}

	/**
	 * This method allows to get the value of the initial velocity passed in the
	 * constructor. If not value is given the default value is zero.
	 * 
	 * @return
	 */

	public double getInitialVelocity() {
		return initialVelocity;
	}

	/**
	 * This method allows to set the value of the initial velocity.
	 * 
	 * @param initialVelocity
	 *            The initial velocity of the particles
	 */

	public void setInitialVelocity(double initialVelocity) {
		this.initialVelocity = initialVelocity;
	}

	/**
	 * This method allows to obtain the value of the inertia weigh. The inertia
	 * can be interpreted as the fluidity of the medium in which a particle
	 * moves.
	 * 
	 * @return The value of the inertia weigh
	 */

	public double getW() {
		return w;
	}

	/**
	 * This method allows to set the inertia weigh value. The inertia can be
	 * interpreted as the fluidity of the medium in which a particle moves.
	 * 
	 * @param w
	 *            The value of the inertia weight
	 */

	public void setW(double w) {
		this.w = w;
	}
	
	/**
	 * This method returns the lower fitness found in each iteration. Thus, allowing to 
	 * keep track of the performance of the algorithm. 
	 * @return
	 * An array that contains the minimum fitness (the lower cost multicast tree) found in each iteration
	 */

	public ArrayList<Double> getMinFitness() {
		return minFitness;
	}

		
	/**
	 * The method allows to add the routing table to the array called
	 * routingTablesArray. It will hold all the routing tables to the destination nodes.
	 * The routing table is obtained from the RoutingTable class.
	 * @param routingTable
	 * An ArrayList of Arraylist that contains the routes to one destination node. 
	 * The routing table is built up using the JgraphT k-shortest path algorithm.
	 * http://jgrapht.org/javadoc/org/jgrapht/alg/KShortestPaths.html
	 */
	public void addRoutingTable(ArrayList<ArrayList<Integer>> routingTable) {
		routingTablesArray.add(routingTable);
	}
	
	/**
	 * The method allows to add the cost matrix. It is obtained from the class GraphCreator
	 * which read a file with the matrix.
	 * @param costMatrix
	 * An array with the cost matrix
	 */

	public void addCostMatrix(int[][] costMatrix) {
		i_costMatrix = new int[costMatrix.length][costMatrix[0].length];
		for (int i = 0; i < costMatrix.length; i++) {
			for (int j = 0; j < costMatrix[i].length; j++) {
				i_costMatrix[i][j] = costMatrix[i][j];
			}
		}
	}
	
	/**
	 * The method allows to add the delay matrix. It is obtained from the class GraphCreator
	 * which read a file with the matrix.
	 * @param delayMatrix
	 * An array with the delay matrix
	 */

	public void addDelayMatrix(int[][] delayMatrix) {
		i_delayMatrix = new int[delayMatrix.length][delayMatrix[0].length];
		for (int i = 0; i < delayMatrix.length; i++) {
			for (int j = 0; j < delayMatrix[i].length; j++) {
				i_delayMatrix[i][j] = delayMatrix[i][j];
			}
		}
	}

	// Here the PSO heuristic starts....

	/**
	 * This method start the PSO heuristic. It generates the initial swarm
	 * randomly. The size of the array that stores the routing tables is taken
	 * as the particle length. That is, the number of destination nodes. Before
	 * using this method is necessary to call the addRoutingTable method times
	 * the number of destination nodes in the multicast communication. In order
	 * to generate the initial swarm the size of the routing tables is taken and
	 * random numbers among this size are generated which represents the index
	 * of the routing table's entries.
	 */

	public void generateInitialSwarm() {
		particleLength = routingTablesArray.size();
		for (int k = 0; k < i_swarmSize; k++) {
			swarm.add(new ArrayList<Integer>());
			pbest.add(new ArrayList<Double>());// Initialized ArrayList of
												// ArrayList to the number of
												// particles
			velocity.add(initialVelocity);
			for (int i = 0; i < particleLength; i++) {
				int routeLength = routingTablesArray.get(i).size();
				int randomRoute = (int) (Math.random() * (routeLength) + 1); // Generate
																				// a
																				// random
																				// number
																				// between
																				// 1
																				// and
																				// the
																				// route
																				// length
				randomRoute = randomRoute - 1;
				swarm.get(k).add(randomRoute);
			}
		}
	}

	/**
	 * Method: calculateVelocity The method calculates the velocity based on the
	 * equation v(t+1) = (w*v)+(c1*rand1*(pbest-x(t)))+(c2*rand2*(gbest-x(t))).
	 * It also checks if the new velocity calculated is within the allowed
	 * bounds. In this case our lower bound is limited to zero (0). This is
	 * because the lower cost of a route cannot be less than zero
	 */

	private void calculateVelocity() {
		rand1 = Math.random(); // never 1
		rand2 = Math.random();
		// select the best position ever by the particle i
		for (int i = 0; i < i_swarmSize; i++) {
			double pbest_i = pbest.get(i).get(0);
			for (int j = 0; j < pbest.get(i).size(); j++) {
				if (pbest_i > pbest.get(i).get(j)) {
					pbest_i = pbest.get(i).get(j);
				}
			}

			newVelocity = (w * velocity.get(i))
					+ (c1 * rand1 * (pbest_i - swarmFitness.get(i)))
					+ (c2 * rand2 * (gbest - swarmFitness.get(i)));
			if (newVelocity < lowerBound) {
				newVelocity = lowerBound;// ***************
			}
			velocityUpdate.add(newVelocity);
		}
	}

	/**
	 * Method: calculatePosition This method update the particle's position. In
	 * order to update the position the equation x(t+1) = x(t)+v(t) is used. The
	 * routing algorithm used a variation of particle swarm optimization based
	 * on the bynary model of PSO. In this variation an elite population exists
	 * and the Sigmoid Function is used as the pdf (probability density
	 * function) which allow to decide whether or not a position i of the
	 * particle changes. The process is as follows: A random number is
	 * generated. Then, the velocity of the particle i is taken an passed to the
	 * sigmoid function which is defined as 1/(1+e^(v)). A result from the
	 * function is obtained and compared to the random number. If the random
	 * number is less than the sigmoid's result the position i of the particle
	 * will be changed for the position i of a random selected particle within
	 * the Nb set.
	 */

	private void calculatePosition() {
		for (int i = 0; i < i_swarmSize; i++) {
			double s = calculateSigmoidFunction(velocityUpdate.get(i));
			for (int j = 0; j < swarm.get(i).size(); j++) {
				double rho = Math.random();
				if (rho < s) {
					int randomIndex = (int) (Math.random() * (bestParticles) + 1); // Generate
																					// a
																					// random
																					// number
																					// between
																					// 1
																					// and
																					// the
																					// bestParticles
																					// number
					randomIndex = randomIndex - 1;
					swarm.get(i).set(j, swarmBest.get(randomIndex).get(j));
				}
			}
		}
		swarmFitness.clear();
		swarmBest.clear();
		swarmFitnessBest.clear();
		velocity.clear();
		velocity.addAll(velocityUpdate);
		velocityUpdate.clear();
	}

	/**
	 * Method: selectBestParticles This method is used to implement the elite
	 * population of particles. The algorithm takes the Nb (number of best)
	 * particles in order to change the position i of a particle is a random
	 * number is less than the output of the sigmoid's function. In order to
	 * select the best Nb particles the array that contains the fitness of the
	 * particles in the swarm is copied to a temporal array. After that, the
	 * method find the position of the lower fitness (it represents the lower
	 * cost multicast tree). Once the position is found, the particle which
	 * belong to the fitness found is copied to the swarmBest array and the
	 * fitness is set to infinite. Thus, the next time the algorithm tries to
	 * find the lower fitness it will not be repeated. The process is repeated
	 * as many times as the Nb variable.
	 */

	private void selectBestParticles() {
		swarmFitnessBest.addAll(swarmFitness);// It’s necessary to create a new
												// ArrayList because original
												// will be changed (change to
												// 0.0)
		for (int i = 0; i < bestParticles; i++) {
			double lower = swarmFitnessBest.get(0);
			int position = 0;
			for (int j = 0; j < swarmFitnessBest.size(); j++) {
				if (swarmFitnessBest.get(j) < lower) {
					lower = swarmFitnessBest.get(j);
					position = j;
				}
			}
			swarmFitnessBest.set(position, Double.POSITIVE_INFINITY);
			orderedFitnessArray.add(lower);// ******check maybe it is not
											// necessary**********
			swarmBest.add(new ArrayList<Integer>());
			for (int j = 0; j < particleLength; j++) {
				swarmBest.get(i).add(swarm.get(position).get(j));
			}

		}
	}

	/**
	 * Method: orderPopulation This method creates a new array called
	 * populationOrderedByFitness. This arrays will contain the the whole swarm
	 * but ordered by its fitness. This is done to check the QoS constrains. The
	 * algorithm tries to find the lower cost multicast tree. When it gets the
	 * last iteration the swarm obtained so far is taken and the QoS constrains
	 * are validated against the particles within. If neither of them satisfies
	 * the QoS constrains then the algorithm reports not to find a multicast
	 * tree with the given QoS constrains. This is done in this way because the
	 * shortest multicast tree has a higher probability to get a lower
	 * end-to-end delay and delay variation-
	 */

	private void orderPopulation() {
		// The whole population ordered according to the fitness. This allow to
		// check QoS constrains the chromosomes
		fitnessTmp.clear();
		orderedFitnessArrayQoS.clear();
		populationOrderedByFitness.clear();
		fitnessTmp.addAll(swarmFitness);// It’s necessary to create a new
										// ArrayList because original will be
										// changed (change to 0.0)
		for (int i = 0; i < i_swarmSize; i++) {
			double lower = fitnessTmp.get(0);
			int position = 0;
			for (int j = 0; j < fitnessTmp.size(); j++) {
				if (fitnessTmp.get(j) < lower) {
					lower = fitnessTmp.get(j);
					position = j;
				}
			}
			fitnessTmp.set(position, Double.POSITIVE_INFINITY);
			orderedFitnessArrayQoS.add(lower);
			populationOrderedByFitness.add(new ArrayList<Integer>());
			for (int j = 0; j < particleLength; j++) {
				populationOrderedByFitness.get(i).add(
						swarm.get(position).get(j));

			}

		}

	}

	/**
	 * Method: calculateLowerFitness This method is used to find the lower
	 * fitness found so far by the whole swarm. This value is used as the gbest
	 * (global best) according to the PSO heuristic. The value gbest is
	 * necessary in order to calculate the new velocity which allows to update
	 * the position of the particles. It also copies the particle holding the
	 * best fitness to the bestParticle array. It is used at the final iteration
	 * to print the lower multicast tree found not necessary that satisfies the
	 * QoS constrains.
	 */

	private void calculateLowerFitness() {
		// calculate fitness
		int position = 0;
		for (int i = 0; i < i_swarmSize; i++) {
			double fitness = calculateFitness(swarm.get(i));
			swarmFitness.add(fitness);
			// System.out.println(fitness);
		}
		for (int i = 0; i < i_swarmSize; i++) {
			pbest.get(i).add(swarmFitness.get(i));
		}

		// find lower fitness
		double lower = swarmFitness.get(0);
		for (int i = 0; i < swarmFitness.size(); i++) {
			if (swarmFitness.get(i) < lower) {
				lower = swarmFitness.get(i);
				position = i;
			}
		}
		gbest = lower;
		bestParticle.clear();
		bestParticle.addAll(swarm.get(position));
		// System.out.println("********************GBest: " + gbest +
		// "*******************");
		minFitness.add(lower);
	}

	/**
	 * Method: calculateFitness This method takes the particle as the input and
	 * extract each position of it. This number is an index of the routing table
	 * corresponding to one of the multicast destination nodes. Once the index
	 * is extracted from the particle, the method takes the array that contains
	 * the routing tables and extract the route (corresponding to the index
	 * extracted). The route is compared with the cost matrix. This allows to
	 * obtain the cost of each of the links that forms the path from source node
	 * to the destination node. The costs are added to obtain the cost of the
	 * path and finally the cost of the multicast tree.
	 * 
	 * @param particle
	 *            A particle of the swarm. It represents a possible solution of
	 *            the problem, that is a multicast routing tree.
	 * @return The fitness of the multicast tree. That is the cost the tree.
	 */

	private double calculateFitness(ArrayList<Integer> particle) {
		double fitness = 0;
		// double sum = 0;
		for (int i = 0; i < particleLength; i++) {
			int counter = 0;
			int sourceNode = 0;
			int destinationNode = 0;
			int route = particle.get(i);
			// String path = routingTablesArray.get(i).get(route).toString();
			// System.out.println(path);
			for (int j = 0; j < routingTablesArray.get(i).get(route).size(); j++) {
				if (j % 2 == 0) {
					sourceNode = routingTablesArray.get(i).get(route).get(j);
					counter++;
				}

				if (j % 2 == 1) {
					destinationNode = routingTablesArray.get(i).get(route)
							.get(j);
					counter++;
				}

				if (counter % 2 == 0) {
					fitness += i_costMatrix[sourceNode][destinationNode];

				}
			}
		}
		return fitness;
	}

	/**
	 * Method: checkDemand This method takes each one of the particles making up
	 * the swarm and checks if one of them satisfies the QoS constrains. In the
	 * routing algorithm this is end-to-end delay and delay variation.
	 * 
	 * @param particle
	 *            The particle containing the indexes of the routes.
	 * @return True if the particle satisfies the QoS constrains false else.
	 */

	private boolean checkDemand(ArrayList<Integer> particle) {
		double delay = 0;
		double delaymi = 0;
		double delaymj = 0;
		double jitter;
		double delaytmp = 0;
		// double sum = 0;
		for (int i = 0; i < particle.size(); i++) {
			int counter = 0;
			int sourceNode = 0;
			int destinationNode = 0;
			int route = particle.get(i);
			// String path = routingTablesArray.get(i).get(route).toString();
			// System.out.println(path);

			for (int j = 0; j < routingTablesArray.get(i).get(route).size(); j++) {
				if (j % 2 == 0) {
					sourceNode = routingTablesArray.get(i).get(route).get(j);
					counter++;
				}

				if (j % 2 == 1) {
					destinationNode = routingTablesArray.get(i).get(route)
							.get(j);
					counter++;
				}

				if (counter % 2 == 0) {
					delay += i_delayMatrix[sourceNode][destinationNode];
					delaytmp += i_delayMatrix[sourceNode][destinationNode];
				}
			}

			if (i == 0) {
				delaymi = delaytmp;
			}
			if (i == 1) {
				delaymj = delaytmp;
			}
			delaytmp = 0;
		}
		jitter = delaymi - delaymj;
		if (jitter <= delayVariation && delay <= delayMulticast) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Method: calculateSigmoidFunction This method takes the velocity of the
	 * particle as the input and based on it the sigmoid's function is used to
	 * calculate a value which will be used as the probability to change the
	 * position i of a particle.
	 * 
	 * @param velocity
	 *            The velocity of the particle.
	 * @return The result of calculation sigmoid's function using the velocity
	 *         as parameter.
	 */

	private double calculateSigmoidFunction(double velocity) {
		double sigmoid = (1 / (1 + (Math.exp(-velocity))));
		return sigmoid;
	}

	/**
	 * Method: getRouteFromParticle This method fills a string with the
	 * multicast tree that satisfies the QoS constrains
	 * 
	 * @param particle
	 *            The particle of the swarm. It represents a possible solution
	 *            of the problem.
	 */

	private void getRouteFromParticle(ArrayList<Integer> particle) {
		for (int i = 0; i < particle.size(); i++) {
			int route = particle.get(i);
			lowerPathQoS += routingTablesArray.get(i).get(route).toString()
					+ " ";
		}
	}

	/**
	 * Method: getLowerCostPath This method fills a string with the lower cost
	 * multicast tree not necessary that satisfies the QoS constrains
	 * 
	 * @param particle
	 *            The particle of the swarm. It represents a possible solution
	 *            of the problem.
	 */

	private void getLowerCostPath(ArrayList<Integer> particle) {
		for (int i = 0; i < particle.size(); i++) {
			int route = particle.get(i);
			lowerPath += routingTablesArray.get(i).get(route).toString() + " ";
		}
	}

	/**
	 * This method is the general loop which performs the iteration of the PSO
	 * heuristic applied to the routing algorithm. When the last iteration is
	 * reached the method calls the checkDemand method in order to check if a
	 * multicast tree that satisfies the QoS constrains was found.
	 */

	public void performParticleRouting() {
		for (int i = 0; i < i_iterations; i++) {
			calculateLowerFitness();
			calculateVelocity();
			selectBestParticles();
			orderPopulation();
			calculatePosition();
//			try {
//				if (i == i_iterations - 1) {
//					boolean qos;
//					int counter = 0;
//					do {
//						qos = checkDemand(populationOrderedByFitness
//								.get(counter));
//						counter++;
//					} while (qos == false);
//					getRouteFromParticle(populationOrderedByFitness
//							.get(counter - 1));
//				}
//
//			} catch (Exception e) {
//				pathFound = false;
//			}
			if (i == i_iterations - 1) {
				getLowerCostPath(bestParticle);
			}

		}
	}

	/**
	 * Method: printLowerCostPathQoS This method shows in the screen the
	 * variable containing the lower multicast tree that satisfies the QoS
	 * constrains found by the routing algorithm.
	 */
	public void printLowerCostPathQoS() {
		if (pathFound == true) {
			System.out.println();
			System.out
					.println("Multicast Tree Found with QoS constrains: (Delay <= "
							+ delayMulticast
							+ " - Delay Variation <= "
							+ delayVariation + ")");
			System.out.println(lowerPathQoS);
		} else {
			System.out
					.println("No Multicast tree found with QoS constrains (Delay <= "
							+ delayMulticast
							+ " - Delay Variation <= "
							+ delayVariation + ")");
		}
	}

	/**
	 * Method: printLowerCostPath This method shows in the screen the variable
	 * containing the lower multicast tree (not necessary satisfying the QoS
	 * constrains) found by the routing algorithm.
	 */
	public void printLowerCostPath() {
		System.out
				.println("The Lower Multicast Tree found by the routing algorithm is: ");
		System.out.println(lowerPath);
	}

	// instance variables
	private double w = 0.729; // inertia weight
	private double c1 = 1.49445; // cognitive weight
	private double c2 = 1.49445; // social weight
	private double initialVelocity = 0;
	private double lowerBound = 0;
	private int bestParticles = 20;
	private int i_swarmSize;
	private int i_iterations;
	private ArrayList<ArrayList<ArrayList<Integer>>> routingTablesArray = new ArrayList<ArrayList<ArrayList<Integer>>>();
	private int[][] i_costMatrix;
	private int[][] i_delayMatrix;
	private int particleLength;
	private ArrayList<ArrayList<Integer>> swarm = new ArrayList<ArrayList<Integer>>();
	private ArrayList<ArrayList<Integer>> swarmBest = new ArrayList<ArrayList<Integer>>();
	private ArrayList<Integer> bestParticle = new ArrayList<Integer>();
	private ArrayList<Double> swarmFitness = new ArrayList<Double>();
	private ArrayList<Double> swarmFitnessBest = new ArrayList<Double>();
	private ArrayList<Double> minFitness = new ArrayList<Double>();
	private double gbest;
	private double rand1; // random values
	private double rand2; // random values
	private ArrayList<ArrayList<Double>> pbest = new ArrayList<ArrayList<Double>>();
	private double newVelocity;
	private ArrayList<Double> velocity = new ArrayList<Double>();
	private ArrayList<Double> velocityUpdate = new ArrayList<Double>();
	private ArrayList<Double> orderedFitnessArray = new ArrayList<Double>();
	private ArrayList<Double> fitnessTmp = new ArrayList<Double>();
	private ArrayList<Double> orderedFitnessArrayQoS = new ArrayList<Double>();
	private ArrayList<ArrayList<Integer>> populationOrderedByFitness = new ArrayList<ArrayList<Integer>>();
	private double delayMulticast = 12;
	private double delayVariation = 5;
	private String lowerPathQoS = "";
	private String lowerPath = "";
	private boolean pathFound = true;
}
