
import java.util.ArrayList;
import java.util.Collections;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author fonty
 */
public class Stats {
    public float mean;
    public float variance;
    public float stdev;
    public int mode;
    public ArrayList<Integer> quantiles;
    public int numQuantiles;
    
    public Stats(ArrayList<Integer> data, int numQuantiles)
    {
        quantiles = new ArrayList<>();
        Collections.sort(data);
        mean = 0;
        this.numQuantiles = numQuantiles;
        float dataPerQuantile = (float) data.size()/numQuantiles;
        
        
        int maxRepeat  = 0;
        int currRepeat = 1;
        for(int i=0; i<data.size(); i++)
        {
            //calculate mode
            if( i > 0 )
            {
                if(data.get(i) == data.get(i-1) )
                {
                    currRepeat++;
                }
                else
                {
                    if(maxRepeat < currRepeat)
                    {
                        maxRepeat  = currRepeat;
                        mode = data.get(i-1);
                    }
                    currRepeat = 0;
                }
            }
            
            mean+= (float) data.get(i)/data.size();
            if(i%dataPerQuantile < 1)
            {
                quantiles.add(data.get(i));
            }
        }
        while(quantiles.size() < numQuantiles)
        {
            quantiles.add(data.get(data.size()-1));
        }
        
        if(currRepeat > maxRepeat)
        {
            mode = data.get(data.size()-1);
        }
        variance = 0;
        
        for(int i=0; i<data.size(); i++)
        {
            variance +=  Math.pow(data.get(i) - mean , 2) / ( data.size() - 1 );
        }
        
        stdev = (float) Math.sqrt(variance);
    }
    
    public String printable()
    {
        String s = mean + "," + stdev + "," + mode;
        
        for(int i=0; i< numQuantiles; i++)
        {
            s+= "," + quantiles.get(i);
        }
        
        
        return s;
    }
}