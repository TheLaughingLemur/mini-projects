import java.io.*;
import java.util.*;
import java.text.*;
import java.math.*;
import java.util.regex.*;
import java.awt.Point;

class ForbiddenTour{
    // grid of current topology of forest
    private char[][] grid;
    // wandCount = number of wand uses;
    public int wandCount;
    // Possible moves that can be made by Hermione and Ron
    private static final Point[] MOVESET = new Point[]{ 
        new Point(1,0),
        new Point(-1,0),
        new Point(0,1),
        new Point(0,-1)
    };
    
    //Constructor
    public ForbiddenTour(char[][] grid){
        this.grid = grid;
        wandCount = 0;
    }
    // tourFrom takes row and col you start from and returns the path to the portkey
    public boolean tourFrom(int row, int col){     
        //check if it's the portkey "*"
        if (grid[row][col] == '*'){
            return true;
        }
           // current position 'M' will now be an "X" as in "a tree grows after we leave a place because we will "never" return;
        grid[row][col] = 'X';
        
        // for each move in the moveset, count up how many possible moves can actually be made;  
        // Also remove them from the moveset for this position I will call "possibleMoves";
        
        int possibleMovesCount = 0;
        Set<Point> possibleMoves = new HashSet<Point>();
        for (Point p: MOVESET){
            int nextRow = row + p.x;
            int nextCol = col + p.y;
            if (nextRow < 0 || nextRow >= grid.length){
                continue;
            }
            if (nextCol < 0 || nextCol >= grid[0].length){
                continue;
            }
            if (grid[nextRow][nextCol] == 'X'){
                continue;
            }
            possibleMovesCount++;        
            possibleMoves.add(p);
        }
     
        
        //Now go through each move and see if a path to portkey can be found;
        if (possibleMoves.size() > 1){
            wandCount++;
        }
        for (Point p: possibleMoves){
            int nextRow = row + p.x;
            int nextCol = col + p.y;
            if(tourFrom(nextRow, nextCol)){
                return true;
            }            
        }
        
        if ( possibleMoves.size() > 1){
            wandCount--;
        }
        return false;  
        
    }   
}



public class Solution {

    public static void main(String[] args) {
        Scanner in = new Scanner(System.in);
        int T = in.nextInt();
        for (int i = 0; i < T; i++){
            int startRow = 0;
            int startCol = 0;
            int n = in.nextInt();
            int m = in.nextInt();
            char[][] forest = new char[n][m];
            for (int row = 0; row < n; row++){
                String line = in.next();
                for (int col = 0; col < m; col++){
                    char item = line.charAt(col);
                    forest[row][col] = item;
                    if(item == 'M'){
                        startRow = row;
                        startCol = col;
                    }
                }
            }
            int k = in.nextInt();
            ForbiddenTour tour = new ForbiddenTour(forest);
            tour.tourFrom(startRow, startCol);  
            int spells = tour.wandCount;
            if(k == spells){
                System.out.println("Impressed");
            } else {
                System.out.println("Oops!");
            }                            
            
        } 
    }
}