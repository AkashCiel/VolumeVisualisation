/*
 * Anna Vilanova: Basic class to represent the gradient. 
 * NO MODIFICATION NEEDED 
 *
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package tudelft.cgv.volume;

/**
 *
 * @author michel
 */
public class VoxelGradient {

    // contais gradient x,y,z component and mag. the magnitude. 
    // the magnitude avoids calculations during the rendering.
    public float x, y, z;
    public float mag;
    
    public VoxelGradient() {
        x = y = z = mag = 0.0f;
    }
    
    public VoxelGradient(float gx, float gy, float gz) {
        x = gx;
        y = gy;
        z = gz;
        mag = (float) Math.sqrt(x*x + y*y + z*z);
    }
    
}
