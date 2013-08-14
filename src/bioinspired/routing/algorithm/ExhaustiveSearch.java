package bioinspired.routing.algorithm;

import java.util.ArrayList;

/**
 * This class is an implementation of the exhaustive search approach. It is used
 * to compare the performance of the routing algorithm with PSO again exhaustive
 * search. The exhaustive search approach is to search one by one all the
 * possible solutions and evaluate them in order to find the lower cost.
 * 
 * @author YEISON JULIAN CAMARGO BARAJAS
 */

public class ExhaustiveSearch {

	/*
	 * Constructor
	 */
	public ExhaustiveSearch() {
	}

	/**
	 * The method allows to add the cost matrix. It is obtained from the class
	 * GraphCreator which read a file with the matrix.
	 * 
	 * @param costMatrix
	 *            An array with the cost matrix
	 */
	public void setCostMatrix(int[][] costMatrix) {
		i_costMatrix = new int[costMatrix.length][costMatrix[0].length];
		for (int i = 0; i < costMatrix.length; i++) {
			for (int j = 0; j < costMatrix[i].length; j++) {
				i_costMatrix[i][j] = costMatrix[i][j];
			}
		}
	}

	/**
	 * The method is the implementation of the exhaustive search.
	 * It receives an ArrayList containing the k-shortest paths to the defined destination nodes.
	 * The k-shortest paths are found with the Routing Table class. Also it receives the destination
	 * nodes and the number of k-shortest paths. 
	 * 
	 * @param routingTable
	 *            An ArrayList with the k-shortest path to the given destination paths.
	 *            The routingTable is obtained from the RoutingTable class with the
	 *            getRoutingTable() method.
	 */
	public double findLowerCostPath(ArrayList<ArrayList<Integer>> routingTable) {

		double lowerCost = Double.POSITIVE_INFINITY;
		int position = 0;
		for (int i = 0; i < routingTable.size(); i++) {
			double cost = 0;
			int counter = 0;
			int sourceNode = 0;
			int destinationNode = 0;
			for (int j = 0; j < routingTable.get(i).size(); j++) {	
				// System.out.println(routingTable.get(i).get(j));

				if (j % 2 == 0) {
					sourceNode = routingTable.get(i).get(j);
					counter++;
				}

				if (j % 2 == 1) {
					destinationNode = routingTable.get(i).get(j);
					counter++;
				}

				if (counter % 2 == 0) {
					cost += i_costMatrix[sourceNode][destinationNode];

				}

			}

			if (cost < lowerCost) {
				lowerCost = cost;
				position = i;

			}
		}

		// System.out.println("Cost: " + lowerCost + "Path: " +
		// routingTable.get(position).toString());

		return 0;

	}

	// instance variables
	private int[][] i_costMatrix;

}
