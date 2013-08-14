package bioinspired.routing.run;

import bioinspired.routing.algorithm.ExhaustiveSearch;
import bioinspired.routing.graph.GraphCreator;
import bioinspired.routing.routingTable.RoutingTable;

public class ExhaustiveSearchRun {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		for (int i = 0; i < 100; i++) {
			
		
		long start = System.currentTimeMillis();
		
		ExhaustiveSearch exhaustive = new ExhaustiveSearch();
		
		GraphCreator g = new GraphCreator("data/topology50NodesFull.txt", 
				"data/CostMatrix50Nodes(1-20).txt", "data/delayMatrix50Nodes(1-5).txt");
		exhaustive.setCostMatrix(g.getCostMatrix());
		
		RoutingTable rTable = new RoutingTable(g.getGraph());
		
		rTable.findKShortestPaths(0, 17, 200);
		exhaustive.findLowerCostPath(rTable.getRoutingTable());
		
        rTable.findKShortestPaths(0, 18, 200);
        exhaustive.findLowerCostPath(rTable.getRoutingTable());
       
        rTable.findKShortestPaths(0, 19, 200);
        exhaustive.findLowerCostPath(rTable.getRoutingTable()); 
        
        long end = System.currentTimeMillis();
        
       
        //double currentMemory = ( (double)((double)(Runtime.getRuntime().totalMemory()/1024)/1024))- ((double)((double)(Runtime.getRuntime().freeMemory()/1024)/1024));
        System.out.println(end-start);

        
		}
  
	}

}
