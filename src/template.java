/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Fonty
 */
public class template {
    String sequence;
    int[] rcc;
    int start;
    int end;
    public template()
    {
        rcc = new int[26];
    }
    public String getRcc()
    {
        String s = ""+rcc[0];
        for(int i=1;i<26;i++)
        {
            s += "," + rcc[i];
        }
        return s;
    }
}
