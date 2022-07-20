import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.junit.Test;

import java.util.ArrayList;

import static org.junit.Assert.*;

public class CrystalSumTest {
    CrystalSum cs = new CrystalSum();

    @Test
    public void dipolarTest() {
        RealVector ri = new ArrayRealVector(new double[]{cs.a, cs.a, cs.c});
        RealVector rj = new ArrayRealVector(new double[]{-cs.a, -2*cs.a, 0});
        assertEquals("Should be equal ", -0.00999041464441944, cs.dipolar(0,2,ri,rj), 1.0e-6);
        assertEquals("Should be equal ", 0.00350255639204227, cs.dipolar(2,2,ri,rj), 1.0e-6);
    }

    @Test
    public void getLocationTest() {
        // verify all spins are counted
        int Lx=5, Ly=8, Lz=2;
        ArrayList<RealVector> spins = new ArrayList<>();
        for (int i=0; i < cs.num_in_cell*8*Lx*Ly*Lz; i++) {
            spins.add(cs.getLocation(i, new int[]{Lx, Ly}));
        }

        for (int i=-Lx; i < Lx; i++){
            for (int j=-Ly; j < Ly; j++){
                for (int k=-Lz; k < Lz; k++){
                    for (int l=0; l < cs.num_in_cell; l++){
                        RealVector r = new ArrayRealVector(3);
                        r.combineToSelf(1, i, cs.latticeVectors[0]);
                        r.combineToSelf(1, j, cs.latticeVectors[1]);
                        r.combineToSelf(1, k, cs.latticeVectors[2]);
                        r.combineToSelf(1, cs.basis[l][0], cs.latticeVectors[0]);
                        r.combineToSelf(1, cs.basis[l][1], cs.latticeVectors[1]);
                        r.combineToSelf(1, cs.basis[l][2], cs.latticeVectors[2]);

                        assertTrue("A spin that should have been was not added by getLocation()", spins.contains(r));
                        spins.remove(r);
                    }
                }
            }
        }
        assertTrue("List should be empty, indicating all spins that were added" +
                "by getLocation() were found directly as well.", spins.isEmpty());

    }

    @Test
    public void directSumTest() {

    }

    @Test
    public void calcPairTest() {
    }

    @Test
    public void calcPairSphereTest() {
    }

}