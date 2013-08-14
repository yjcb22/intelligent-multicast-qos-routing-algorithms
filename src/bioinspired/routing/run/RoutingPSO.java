	package bioinspired.routing.run;
	
	
	
	
	import java.util.Random;
	
	import bioinspired.routing.algorithm.RoutingParticleSwarm;
	import bioinspired.routing.graph.GraphCreator;
	import bioinspired.routing.matrix.Matrix;
	import bioinspired.routing.performance.GeneticMetaOptimizer;
import bioinspired.routing.routingTable.RoutingTable;
	
	public class RoutingPSO {
		
		public static void main (String[] args) {
			
			for (int i = 0; i < 100; i++) {
			long start = System.currentTimeMillis();
			
			GraphCreator g = new GraphCreator("data/topology30NodesFull.txt", 
					"data/CostMatrix30Nodes(1-20).txt", "data/delayMatrix30Nodes(1-5).txt");		
			//RoutingParticleSwarm rp = new RoutingParticleSwarm(40, 15);
			RoutingParticleSwarm rp = new RoutingParticleSwarm(46, 11, 0.5258424935, 1.0548425422, 0.0625287829, 14, 10, 15, 10);
			RoutingTable rTable = new RoutingTable(g.getGraph());
			
			rTable.findKShortestPaths(0, 17, 200);	
	        rp.addRoutingTable(rTable.getRoutingTable());
	        rTable.findKShortestPaths(0, 18, 200);
	        rp.addRoutingTable(rTable.getRoutingTable());
	        rTable.findKShortestPaths(0, 19, 200);
	        rp.addRoutingTable(rTable.getRoutingTable());

			rp.addCostMatrix(g.getCostMatrix());
			//rp.addDelayMatrix(g.getDelayMatrix());

			rp.generateInitialSwarm();
			rp.performParticleRouting();
			
			long end = System.currentTimeMillis();		        
	        System.out.println(end-start);	
	        
//	        rp.printLowerCostPath();
//	        rp.printLowerCostPathQoS();
//	        System.out.println(rp.getMinFitness().toString());
				
				
			//double currentMemory = ( (double)((double)(Runtime.getRuntime().totalMemory()/1024)/1024))- ((double)((double)(Runtime.getRuntime().freeMemory()/1024)/1024));
		       

		        

				
			}			
			
			
		}
	
	}
