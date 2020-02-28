/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Fonty
 */
public class Domain{
    int numSegments;
    String[] chainStart;
    String[] chainEnd;
    int[] start;
    int[] end;
    String name;
    Domain(int numSegments)
    {
        this.numSegments = numSegments;
        this.chainStart = new String[numSegments];
        this.chainEnd = new String[numSegments];
        this.start = new int[numSegments];
        this.end   = new int[numSegments];
        name = "";
    }
}