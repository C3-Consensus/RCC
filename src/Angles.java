
import java.util.ArrayList;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author fonty
 */
public class Angles {
    public double psi;
    public double phi;
    public double omega;
    
    public Angles()
    {
        this.phi = 0;
        this.psi = 0;
        this.omega = 0;
    }
    
    public Angles(double psi, double phi, double omega)
    {
        this.psi = psi;
        this.phi = phi;
        this.omega = omega;
    }
    
    public Angles(Angles t)
    {
        this.psi = t.psi;
        this.phi = t.phi;
        this.omega = t.omega;
    }
    
    public static double samplePrior(ArrayList<Double> prior)
    {
        int quantile = (int) Math.round(Constants.rand.nextDouble()*1000);
        return prior.get(quantile);
    }
    
    public static AngleSample sampleLocalPrior(ArrayList<Double> prior, double s0, int numNeighbours)
    {
        int quantile = prior.size()/2;
        //find belonging quantile for s0
        int searchStep = quantile;
        for(int i=0; i< Constants.numSearches; i++)
        {
            searchStep /= 2;
//        System.out.println(quantile);
            if(prior.get(quantile) == s0)
            {
                break;
            }
            if( prior.get(quantile) > s0)
            {
                quantile -= searchStep;
            }
            else
            {
                quantile += searchStep;
            }
        }
        
        quantile += Constants.rand.nextInt(numNeighbours) - ( numNeighbours/2 );
        if(quantile < 0)
        {
            quantile += Constants.numQuantiles;
        }
        else if(quantile > Constants.numQuantiles)
        {
            quantile -= Constants.numQuantiles;
        }
        
        return new AngleSample(prior.get(quantile), quantile);
    }
    
    public static AngleSample sampleLocalPrior(ArrayList<Double> prior, AngleSample s0, int numNeighbours)
    {
        int quantile = s0.quantile + Constants.rand.nextInt(numNeighbours) - (numNeighbours/2);
        if(quantile<0)
        {
            quantile += Constants.numQuantiles;
        }
        else if(quantile> Constants.numQuantiles)
        {
            quantile -= Constants.numQuantiles;
        }
        return new AngleSample(prior.get(quantile), quantile);
    }
}
