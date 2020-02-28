
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

public class Reticle{
    public double minX;
    public double minY;
    public double minZ;
    public double maxX;
    public double maxY;
    public double maxZ;
    public double distance;
    public double sideLength;
    public int numCubesX;
    public int numCubesY;
    public int numCubesZ;
    public double maxRadius;

    public ArrayList<ArrayList<Integer>> cubeMembers;

    public Reticle( double minX, double maxX, double minY, double maxY, double minZ, double maxZ)
    {
        //initialize bounds of reticle
        this.minX = minX;
        this.minY = minY;
        this.minZ = minZ;
        this.maxX = maxX;
        this.maxY = maxY;
        this.maxZ = maxZ;
    }

    Reticle(Reticle r)
    {
        minX = r.minX;
        minY = r.minY;
        minZ = r.minZ;
        maxX = r.maxX;
        maxY = r.maxY;
        maxZ = r.maxZ;
        distance = r.distance;
        sideLength = r.sideLength;
        numCubesX = r.numCubesX;
        numCubesY = r.numCubesY;
        numCubesZ = r.numCubesZ;
        maxRadius = r.maxRadius;    
    }

    public void Initialize(double distance)
    {
        //fix side length for each cube
        //this.distance = distance + 3* maxRadius;
        this.sideLength = distance + 2 * maxRadius + .001; //this.distance / Math.sqrt(3);

        //calculate number of cubes per side
        this.numCubesX = (int) (Math.floor((maxX - minX )/sideLength ) + 1) ;
        this.numCubesY = (int) (Math.floor((maxY - minY )/sideLength ) + 1) ;
        this.numCubesZ = (int) (Math.floor((maxZ - minZ )/sideLength ) + 1) ;

        //initialize the list of contents for each cube with nulls  
        cubeMembers = new ArrayList<>();
        for(int i = 0 ; i < numCubesX * numCubesY * numCubesZ ; i++)
        {
            cubeMembers.add(null);
        }

    }


    public Reticle( double minX, double maxX, double minY, double maxY, double minZ, double maxZ, double distance)
    {
        //initialize bounds of reticle
        this.minX = minX;
        this.minY = minY;
        this.minZ = minZ;
        this.maxX = maxX;
        this.maxY = maxY;
        this.maxZ = maxZ;

        //fix side length for each cube
        this.distance = distance;
        this.sideLength = distance / Math.sqrt(3);

        //calculate number of cubes per side
        this.numCubesX = (int) (Math.floor((maxX - minX )/sideLength ) + 1) ;
        this.numCubesY = (int) (Math.floor((maxY - minY )/sideLength ) + 1) ;
        this.numCubesZ = (int) (Math.floor((maxZ - minZ )/sideLength ) + 1) ;

        //initialize the list of contents for each cube with nulls  
        cubeMembers = new ArrayList<>();
        for(int i = 0 ; i < numCubesX * numCubesY * numCubesZ ; i++)
        {
            cubeMembers.add(null);
        }
    }

    public int HashCube(double x, double y, double z)
    {
        return (int) (Math.floor( (x - minX) / this.sideLength) 
                    + Math.floor( (y - minY) / this.sideLength) * this.numCubesX
                    + Math.floor( (z - minZ) / this.sideLength) * this.numCubesX * this.numCubesY);
    }

    //return the residues id in neighbour cubes excluding the current one
    public ArrayList<Integer> getNeighbours (int c)
    {
        ArrayList<Integer> neighbours = new ArrayList<>();
        ArrayList<Integer> dummy = new ArrayList<>();

        int nxy = numCubesX * numCubesY;
        int nxy1 = numCubesX * (numCubesY -1);
        int nxyz1 = numCubesX * numCubesY * (numCubesZ-1);

        for( int x = -1 ; x < 2; x++)
        {
        //if cube is in left border in x axis
            if( x == -1 )
            {
                if( c % numCubesX == 0)
                {
                    continue;
                }
            }
        //if cube is in right border in x axis                    
            if( x == 1)
            {
                if( c % numCubesX == numCubesX-1)
                {
                    continue;
                }
            }

            for( int y = -1; y < 2; y++)
            {
            //if cube is in left border in y axis
                if( y == -1 )
                {
                    if( c % nxy < numCubesX )
                    {
                        continue;
                    }
                }
            //if cube is in right border in y axis
                if( y == 1)
                {
                    if( c % nxy >= nxy1)
                    {
                        continue;
                    }
                }

                for( int z = -1; z < 2; z++)
                {
                //if cube is in left border in z axis
                    if(z == -1 )
                    {
                        if( c < nxy)
                        {
                            continue;
                        }
                    }
                //if cube is in left border in z axis
                    if( z == 1)
                    {
                        if( c >= nxyz1)
                        {
                            continue;
                        }
                    }
                    //              c +- 1 +-   numCubesX    +- numCubesX * numCubesY
                    int neighbour = c +  x +   y* numCubesX  + z * numCubesX * numCubesY;
                    dummy.add(neighbour);
                    /*
                    //skip current cube
                    if (neighbour == c)
                    {
                        continue;
                    }
                    */
                    //check if neighbour cube is initialized
                    if( cubeMembers.get(neighbour) != null)
                    {
                        //add residues from current cube to neighbours list
                        for( int r = 0; r< cubeMembers.get(neighbour).size(); r++)
                        {
                            neighbours.add( cubeMembers.get(neighbour).get(r));
                        }
                    }
                }
            }
        }

        return neighbours;
    }

                //return the residues id in neighbour cubes excluding the current one
    public ArrayList<Integer> getNeighboursUpper (int c)
    {
        ArrayList<Integer> neighbours = new ArrayList<>();
        ArrayList<Integer> dummy = new ArrayList<>();

        int nxy = numCubesX * numCubesY;
        int nxy1 = numCubesX * (numCubesY -1);
        int nxyz1 = numCubesX * numCubesY * (numCubesZ-1);

        for( int x = 0 ; x < 2; x++)
        {
        //if cube is in right border in x axis                    
            if( x == 1)
            {
                if( c % numCubesX == numCubesX-1)
                {
                    continue;
                }
            }

            for( int y = 0; y < 2; y++)
            {
            //if cube is in right border in y axis
                if( y == 1)
                {
                    if( c % nxy >= nxy1)
                    {
                        continue;
                    }
                }

                for( int z = 0; z < 2; z++)
                {
                //if cube is in left border in z axis
                    if( z == 1)
                    {
                        if( c >= nxyz1)
                        {
                            continue;
                        }
                    }
                    //              c +- 1 +-   numCubesX    +- numCubesX * numCubesY
                    int neighbour = c +  x +   y* numCubesX  + z * numCubesX * numCubesY;
                    dummy.add(neighbour);
                    /*
                    //skip current cube
                    if (neighbour == c)
                    {
                        continue;
                    }
                    */
                    //check if neighbour cube is initialized
                    if( cubeMembers.get(neighbour) != null)
                    {
                        //add residues from current cube to neighbours list
                        for( int r = 0; r< cubeMembers.get(neighbour).size(); r++)
                        {
                            neighbours.add( cubeMembers.get(neighbour).get(r));
                        }
                    }
                }
            }
        }

        return neighbours;
    }

}
        
