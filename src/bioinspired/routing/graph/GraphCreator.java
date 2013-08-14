package bioinspired.routing.graph;

import java.io.*;
import java.util.*;

import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.*;
import org.jgrapht.graph.DefaultDirectedGraph;

/**
 * This class is the base for the routing algorithms. It creates the graph from
 * a text file with adjacency matrix representation. In order to create the
 * graph, the JGraphT (http://jgrapht.org/) library is used.
 * 
 * @author YEISON JULIAN CAMARGO BARAJAS
 * 
 */
public class GraphCreator {

	/**
	 * The constructor of the Graph class. It creates a directed graph from a
	 * file with and adjacency matrix representation of it. This must be the
	 * first class called when performing bio-inspired routing. Based on the
	 * graph created the routing algorithm will create the routing tables and
	 * execute the routing.
	 * 
	 * @param pathToTopologyFile
	 *            The name of the path in the SO to the file to be loaded. It
	 *            has to be an absolute path. The file must be a text file which
	 *            contains a matrix to represent the graph (in the adjacency
	 *            matrix form).
	 * 
	 * @param pathToCostMatrix
	 *            The name of the path in the SO to the file to be loaded. It
	 *            has to be an absolute path. The file must be a text file which
	 *            contains the cost of the link connecting the nodes in the
	 *            topology. The cost matrix is the same of a topology matrix
	 *            (adjacency representation matrix) but instead of having the
	 *            connections between the nodes it contains a number with is the
	 *            cost of the link
	 * 
	 * @param pathToDelayMatrix
	 *            The name of the path in the SO to the file to be loaded. It
	 *            has to be an absolute path. The file must be a text file which
	 *            contains the delay of the links connecting the nodes in the
	 *            topology. The delay matrix is the same of a topology matrix
	 *            (adjacency representation matrix) but instead of having the
	 *            connections between the nodes it contains a number with is the
	 *            delay of the link
	 */
	// Constructor
	public GraphCreator(String pathToTopologyFile, String pathToCostMatrix,
			String pathToDelayMatrix) {
		this.pathToTopologyFile = pathToTopologyFile;
		this.pathToCostMatrix = pathToCostMatrix;
		this.pathToDelayMatrix = pathToDelayMatrix;
		generateGraphFromFile();
		generateCostMatrixFromFile();
		generateDelayMatrixFromFile();
	}

	// Getters
	/**
	 * This method returns the graph created from the topology text file. It
	 * will be used by the rest of the class in the routing algorithm
	 * 
	 * @return The method returns an object of type graph in the jgrapht
	 *         project.
	 */
	public Graph<Long, DefaultEdge> getGraph() {
		return g;
	}

	/**
	 * This method return the adjacency matrix created from the topology text
	 * file. It can be passed to other classes.
	 * 
	 * @return An int matrix which indicates the connection between the nodes of
	 *         the graph
	 */

	public int[][] getAdjacencyMatrix() {
		return adjacencyMatrix;
	}

	/**
	 * This method return the cost matrix created from the cost text file. It
	 * can be passed to other classes.
	 * 
	 * @return An int matrix which indicates the cost of the link connecting the
	 *         nodes of the graph
	 */

	public int[][] getCostMatrix() {
		return costMatrix;
	}

	/**
	 * This method return the delay matrix created from the delay text file. It
	 * can be passed to other classes.
	 * 
	 * @return An int matrix which indicates the delay of the link connecting
	 *         the nodes of the graph
	 */

	public int[][] getDelayMatrix() {
		return delayMatrix;
	}

	/**
	 * Method: generateGraphFromFile This method takes the file where the matrix
	 * representing the topology resides and convert it to a matrix that is used
	 * to create the graph using the JgraphT project
	 */
	private void generateGraphFromFile() {
		int i = 0;
		int j = 0;
		try {
			FileInputStream fstream = new FileInputStream(pathToTopologyFile);
			FileInputStream fstream2 = new FileInputStream(pathToTopologyFile);
			DataInputStream in = new DataInputStream(fstream);
			DataInputStream in2 = new DataInputStream(fstream2);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			BufferedReader br2 = new BufferedReader(new InputStreamReader(in2));
			String strLine;
			String line = br2.readLine();
			StringTokenizer nodeCount = new StringTokenizer(line, " ");
			adjacencyMatrix = new int[nodeCount.countTokens()][nodeCount
					.countTokens()];
			while ((strLine = br.readLine()) != null) {
				// Print the content on the console
				// System.out.println(strLine);
				StringTokenizer nodeId = new StringTokenizer(strLine, " ");
				// System.out.println(nodeId.countTokens());
				int nodes = nodeId.countTokens();

				while (nodeId.hasMoreTokens()) {
					String digits = nodeId.nextToken();
					adjacencyMatrix[i][j] = Integer.parseInt(digits.trim());
					j++;
					// System.out.println(Integer.parseInt(digits.trim()));
					// System.out.println("M: " + adjacencyMatrix[i][j]);
				}
				i++;
				j = 0;
			}
			in.close();
		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
		}
		// vertices = new Object[adjacencyMatrix.length]; ***************

		// System.out.println("Topology Matrix: ");
		// for (int k = 0; k < adjacencyMatrix.length; k++) {
		// for (int l = 0; l < adjacencyMatrix[k].length; l++) {
		// System.out.print(adjacencyMatrix[k][l] + " ");
		// }
		// System.out.println();
		// }

		/**
		 * Once the matrix is created, this is used to create the vertices and
		 * edges of the graph in the JgraphT project.
		 */
		for (int k = 0; k < adjacencyMatrix.length; k++) {
			g.addVertex((long) k);
		}
		for (int k = 0; k < adjacencyMatrix.length; k++) {
			for (int l = 0; l < adjacencyMatrix[k].length; l++) {
				if (adjacencyMatrix[k][l] == 1) {
					g.addEdge((long) k, (long) l);
				}
			}

		}

	}

	// Generate cost matrix, size must be the same of the topology matrix.
	/**
	 * Method: generateCostMatrixFromFile This method takes the file where the
	 * cost of the link connecting the nodes of the graph resides and convert it
	 * to a matrix that is used to be passed to other classes in order to
	 * perform the bio-inspired routing with QoS constrains.
	 */

	private void generateCostMatrixFromFile() {
		int i = 0;
		int j = 0;
		try {
			FileInputStream fstream = new FileInputStream(pathToCostMatrix);
			FileInputStream fstream2 = new FileInputStream(pathToCostMatrix);
			DataInputStream in = new DataInputStream(fstream);
			DataInputStream in2 = new DataInputStream(fstream2);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			BufferedReader br2 = new BufferedReader(new InputStreamReader(in2));
			String strLine;
			String line = br2.readLine();
			StringTokenizer nodeCount = new StringTokenizer(line, " ");
			costMatrix = new int[nodeCount.countTokens()][nodeCount
					.countTokens()];
			while ((strLine = br.readLine()) != null) {
				StringTokenizer nodeId = new StringTokenizer(strLine, " ");
				int nodes = nodeId.countTokens();
				while (nodeId.hasMoreTokens()) {
					String digits = nodeId.nextToken();
					costMatrix[i][j] = Integer.parseInt(digits.trim());
					j++;
				}
				i++;
				j = 0;
			}
			in.close();

		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
		}

	}

	// Generate delay matrix, size must be the same of the topology matrix.
	/**
	 * Method: generateDelayMatrixFromFile This method takes the file where the
	 * delay of the link connecting the nodes of the graph resides and convert
	 * it to a matrix that is used to be passed to other classes in order to
	 * perform the bio-inspired routing with QoS constrains.
	 */

	private void generateDelayMatrixFromFile() {
		int i = 0;
		int j = 0;
		try {
			FileInputStream fstream = new FileInputStream(pathToDelayMatrix);
			FileInputStream fstream2 = new FileInputStream(pathToDelayMatrix);
			DataInputStream in = new DataInputStream(fstream);
			DataInputStream in2 = new DataInputStream(fstream2);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			BufferedReader br2 = new BufferedReader(new InputStreamReader(in2));
			String strLine;
			String line = br2.readLine();
			StringTokenizer nodeCount = new StringTokenizer(line, " ");
			delayMatrix = new int[nodeCount.countTokens()][nodeCount
					.countTokens()];
			while ((strLine = br.readLine()) != null) {
				StringTokenizer nodeId = new StringTokenizer(strLine, " ");
				int nodes = nodeId.countTokens();
				while (nodeId.hasMoreTokens()) {
					String digits = nodeId.nextToken();
					delayMatrix[i][j] = Integer.parseInt(digits.trim());
					j++;
				}
				i++;
				j = 0;
			}
			in.close();

		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
		}
	}

	public void printTopologyMatrix() {
		System.out.println("Topology Matrix: ");
		for (int i = 0; i < adjacencyMatrix.length; i++) {
			for (int j = 0; j < adjacencyMatrix[i].length; j++) {
				System.out.print(adjacencyMatrix[i][j] + " ");
			}
			System.out.println();
		}
	}

	public void printCostMatrix() {
		System.out.println("Cost Matrix: ");
		for (int i = 0; i < costMatrix.length; i++) {
			for (int j = 0; j < costMatrix[i].length; j++) {
				System.out.print(costMatrix[i][j] + " ");
			}
			System.out.println();
		}
	}

	// Instance Variables
	// private String pathToFile;
	private String pathToTopologyFile;
	private String pathToCostMatrix;
	private String pathToDelayMatrix;
	private Graph<Long, DefaultEdge> g = new DefaultDirectedGraph<Long, DefaultEdge>(
			DefaultEdge.class);
	private int[][] adjacencyMatrix;
	private int[][] costMatrix;
	private int[][] delayMatrix;

}
