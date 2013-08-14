package bioinspired.routing.matrix;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.Set;
import java.util.StringTokenizer;
import org.jgrapht.ListenableGraph;
import org.jgrapht.graph.DefaultEdge;

/**
*
* @author yeison
*/
public class Matrix {

   /*
    * Constructor
    */
   public Matrix(org.jgrapht.Graph graph) {
       g = graph;
       NumberOfVertices = g.vertexSet().size();
       edgesAll = g.edgeSet();
       buildTopologyMatrix();
       buildCostMatrix();
   }

   /*
    * Getters
    */
   public int[][] getAdjacencyMatrix() {
       return adjacencyMatrix;
   }

   public int[][] getCostMatrix() {
       return costMatrix;
   }
   
   public int[][] getDelayMatrix() {
       return delayMatrix;
   }

   /*
    * Method : buildTopologyMatrix
    */
   /*
    * Gets a instance of the topology built in the Graph class and generates a
    * representation of the network in the adjacency matrix form.
    */
   private void buildTopologyMatrix() {
       adjacencyMatrix = new int[NumberOfVertices][NumberOfVertices];

       for (DefaultEdge edge : edgesAll) {
           //System.out.println(edge);
           String node = edge.toString();
           int sourceNode = 0;
           int destinationNode = 0;
           boolean destination = false;
           StringTokenizer nodeId = new StringTokenizer(node, ":() ");
           while (nodeId.hasMoreTokens()) {
               String digits = nodeId.nextToken();
               //System.out.println(Integer.parseInt(digits.trim()));
               if (destination == false) {
                   sourceNode = Integer.parseInt(digits.trim());
                   destination = true;
               }
               if (destination == true) {
                   destinationNode = Integer.parseInt(digits.trim());
               }
           }
           adjacencyMatrix[sourceNode][destinationNode] = 1;
       }
   }

   /*
    * Method : buildCostMatrix
    */
   /*
    * Based on the adjacency matrix built with the buildTopologyMatrix method,
    * it generates a random cost matrix.
    */
   private void buildCostMatrix() {
       costMatrix = new int[NumberOfVertices][NumberOfVertices];

       for (int i = 0; i < costMatrix.length; i++) {
           for (int j = 0; j < costMatrix[i].length; j++) {
               int randomCost = (int) (Math.random() * (10) + 1); // Generate Random cost from 1 to 5
               costMatrix[i][j] = randomCost;
           }
       }

       for (int i = 0; i < costMatrix.length; i++) { //Diagonal to zero (loop)
           for (int j = 0; j < costMatrix[i].length; j++) {
               if (i == j) {
                   costMatrix[i][j] = 0;
               }
           }
       }

   }

   /*
    * Method : readCostMatrix
    */
   /*
    * Reads the cost matrix from a file a stores it in a local variable.
    */
   public void readCostMatrix(String pathToFile) {
//http://www.roseindia.net/java/beginners/java-read-file-line-by-line.shtml
       costMatrix = new int[NumberOfVertices][NumberOfVertices];
       int i = 0;
       int j = 0;
       try {

           FileInputStream fstream = new FileInputStream(pathToFile);
           DataInputStream in = new DataInputStream(fstream);
           BufferedReader br = new BufferedReader(new InputStreamReader(in));
           String strLine;

           while ((strLine = br.readLine()) != null) {
               // Print the content on the console
               //System.out.println(strLine);
               StringTokenizer nodeId = new StringTokenizer(strLine, " ");
               while (nodeId.hasMoreTokens()) {
                   String digits = nodeId.nextToken();
                   costMatrix[i][j] = Integer.parseInt(digits.trim());
                   j++;
                   //System.out.println(Integer.parseInt(digits.trim()));
               }
               i++;
               j = 0;
           }
           in.close();

       } catch (Exception e) {
           // TODO: handle exception
           System.err.println("Error: " + e.getMessage());
       }
   }
   
   /*
    * Method : readCostMatrix
    */
   /*
    * Reads the cost matrix from a file a stores it in a local variable.
    */
   public void readDelayMatrix(String pathToFile) {
//http://www.roseindia.net/java/beginners/java-read-file-line-by-line.shtml
       delayMatrix = new int[NumberOfVertices][NumberOfVertices];
       int i = 0;
       int j = 0;
       try {

           FileInputStream fstream = new FileInputStream(pathToFile);
           DataInputStream in = new DataInputStream(fstream);
           BufferedReader br = new BufferedReader(new InputStreamReader(in));
           String strLine;

           while ((strLine = br.readLine()) != null) {
               // Print the content on the console
               //System.out.println(strLine);
               StringTokenizer nodeId = new StringTokenizer(strLine, " ");
               while (nodeId.hasMoreTokens()) {
                   String digits = nodeId.nextToken();
                   delayMatrix[i][j] = Integer.parseInt(digits.trim());
                   j++;
                   //System.out.println(Integer.parseInt(digits.trim()));
               }
               i++;
               j = 0;
           }
           in.close();

       } catch (Exception e) {
           // TODO: handle exception
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
   //instance variables
   private org.jgrapht.Graph g;
   private int NumberOfVertices;
   private Set<DefaultEdge> edgesAll;
   private int[][] adjacencyMatrix;
   private int[][] costMatrix;
   private int[][] delayMatrix;
}