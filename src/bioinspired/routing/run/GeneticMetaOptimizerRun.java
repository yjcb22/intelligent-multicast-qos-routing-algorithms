package bioinspired.routing.run;

import bioinspired.routing.performance.GeneticMetaOptimizer;

public class GeneticMetaOptimizerRun {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		GeneticMetaOptimizer gmo = new GeneticMetaOptimizer(100, 100, 0.1);
		//System.out.println(gmo.getPopulation().toString()+" "+gmo.getPopulation().size());
		System.out.println(gmo.getLowerFitnessPerGeneration().toString());
		System.out.println(gmo.getLowerChromosomePerGeneration().toString());

	}

}
