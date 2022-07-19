import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

public class CrystalSum {


    public static double dipolar(int comp1, int comp2, RealVector ri, RealVector rj){
        double norm = ri.subtract(rj).getNorm();
        if (comp1==comp2){
            return (norm*norm - 3*Math.pow(ri.subtract(rj).getEntry(comp1), 2))/Math.pow(norm, 5);
        } else {
            return (-3*ri.subtract(rj).getEntry(comp1)*ri.subtract(rj).getEntry(comp2))/Math.pow(norm, 5);
        }
    }

    public static RealVector getLocation(int n, int[] limitingBox, int num_in_cell, RealVector[] latticeVectors, double[][] basis){
        int Lx = limitingBox[0], Ly = limitingBox[1];
        int k = ((n/num_in_cell)/Lx)/Ly/8;
        int sym_op = ((n/num_in_cell)/Lx)/Ly%8;
        int j = ((n/num_in_cell)/Lx)%Ly;
        int i = (n/num_in_cell)%Lx;
        int l = n%num_in_cell;

        int[] symVec = new int[]{sym_op%2,sym_op/2%2,sym_op/4%2};

        int[] locCoordinates = new int[]{i, j, k};
        RealVector loc = new ArrayRealVector(3);
        for (int coor=0;coor<latticeVectors.length;coor++){
            loc.combineToSelf(1, locCoordinates[coor], latticeVectors[coor]);   // add primitive lattice vectors according to (i,j,k)
            loc.combineToSelf(1, basis[l][coor], latticeVectors[coor]);         // add fraction of primitive lattice vector according to basis index l
            if (coor != 2 && symVec[coor] != 0) {
                loc.combineToSelf(1, -limitingBox[coor], latticeVectors[coor]);
            }
        }
        if (symVec[2]==1){
            loc = loc.subtract(latticeVectors[2].mapMultiply(2*k+1));
        }
        return loc;
    }

    public static double directSum(RealVector[] latticeVectors, double[][] basis, double heightTol, double radiusTol, int num_in_cell){
        RealVector origin = new ArrayRealVector(3); // constructs a vector of 3 zeros
        double prevSum=0, LzSum=0;
        for (int Lx=1; Lx<300; Lx++){
            int Ly=Lx;
            double prevLzSum=0;
            LzSum=0;
            int i;
            for (i=1;; i++) {
                RealVector ri = getLocation(i, new int[]{Lx, Ly}, num_in_cell, latticeVectors, basis);
                if (ri.getSubVector(0, 2).getNorm() < Lx) {
//                    LzSum += dipolar(2, 2, origin, ri);
                    LzSum += calcPair(origin, ri, 0.00000005, 0.00005, num_in_cell, latticeVectors, basis);
                }
                if (i % (num_in_cell * Lx * Ly * 8) == 0 && i / (num_in_cell * Lx * Ly * 8) > Lx) {    // a new "floor" is added at the top and bottom AND the height of the cylinder is larger than its basis radius
                    System.out.println("Lx\t" + Lx + " Lz\t" + (i / (num_in_cell * Lx * Ly * 8)) + " Lz_sum\t " + LzSum);
                    if (Math.abs((LzSum - prevLzSum)/LzSum) < heightTol) {  // convergence
                        System.out.println("Lz convergence! \t Lx\t" + Lx + " Lz\t" + (i / (num_in_cell * Lx * Ly * 8)) + " Lz_sum\t " + LzSum);
                        break;
                    }
                    prevLzSum = LzSum;
                }
            }
            if (Math.abs((LzSum - prevSum)/LzSum) < radiusTol){
                System.out.println("Lx convergence!! \t Lx\t" + Lx + " Lz\t" + (i / (num_in_cell * Lx * Ly * 8)) + " Lz_sum\t " + LzSum);
                break;
            }
            prevSum=LzSum;
        }

        return LzSum;
    }

    public static double calcPair(RealVector ri, RealVector rj, double heightTol, double radiusTol, int num_in_cell, RealVector[] latticeVectors, double[][] basis){
        double prevSum=0, LzSum;
        for (int Lx = (int)(ri.subtract(rj).getSubVector(0,2).getNorm()+1);;Lx++){
            int Ly = Lx;
            double prevLzSum=0;
            LzSum=0;
            int k;
            for (k=1;;k++){
                RealVector rk = getLocation(k, new int[]{Lx, Ly}, num_in_cell, latticeVectors, basis);
                if (rk.getSubVector(0, 2).getNorm() < Lx && !ri.equals(rk) && !rj.equals(rk)){
                    LzSum += dipolar(0, 2, ri,rk) * dipolar(0, 2, rk,rj);
                }
                if (k % (num_in_cell * Lx * Ly * 8) == 0 && k / (num_in_cell * Lx * Ly * 8) > Lx) {    // a new "floor" is added at the top and bottom AND the height of the cylinder is larger than its basis radius
                    if (Math.abs((LzSum - prevLzSum)/LzSum) < heightTol) {  // convergence
//                        System.out.println("Lz convergence! \t Lx\t" + Lx + " Lz\t" + (k / (num_in_cell * Lx * Ly * 8)) + " Lz_sum\t " + LzSum);
                        break;
                    }
                    prevLzSum = LzSum;
                }
            }
            if (Math.abs((LzSum - prevSum)/LzSum) < radiusTol){
//                System.out.println("Lx convergence!! \t Lx\t" + Lx + " Lz\t" + (k / (num_in_cell * Lx * Ly * 8)) + " Lz_sum\t " + LzSum);
                break;
            }
            prevSum=LzSum;
        }
        return LzSum;
    }

    public static void main(String[] args){
        double a=1;
        double c=2.077294686;
        int num_in_cell=4;
        RealVector[] latticeVectors = new RealVector[3];
        latticeVectors[0] = new ArrayRealVector(new double[]{a, 0, 0});
        latticeVectors[1] = new ArrayRealVector(new double[]{0, a, 0});
        latticeVectors[2] = new ArrayRealVector(new double[]{0, 0, c});
        double[][] basis = {
            {0,0,0},
            {0.5,0,0.25},
            {0.5,0.5,0.5},
            {0,0.5,0.75}
        };

        RealVector origin = new ArrayRealVector(3); // constructs a vector of 3 zeros
        System.out.println(directSum(latticeVectors, basis, 0.005, 0.005, num_in_cell));
    }
}
