/*
 * Anna Vilanova: This function does not need to be modified.  
 * NO MODIFICATION NEEDED 
 *
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package tudelft.cgv.volvis;

import com.jogamp.opengl.GL2;
import java.util.ArrayList;
import tudelft.cgv.util.TFChangeListener;
/**
 *
 * @author michel
 */
public abstract class Renderer {
     int winWidth, winHeight;
    boolean visible = false;
    boolean interactiveMode = false;
    ArrayList<TFChangeListener> listeners = new ArrayList<TFChangeListener>();

    //Do NOT modify this function
    public Renderer() {
        
    }
    //Do NOT modify this function
    public void setInteractiveMode(boolean flag) {
        interactiveMode = flag;
    }
    //Do NOT modify this function
    public void setWinWidth(int w) {
        winWidth = w;
    }
    //Do NOT modify this function
    public void setWinHeight(int h) {
        winHeight = h;
    }
    //Do NOT modify this function
    public int getWinWidth() {
        return winWidth;
    }
    //Do NOT modify this function
    public int getWinHeight() {
        return winHeight;
    }
    //Do NOT modify this function
    public void setVisible(boolean flag) {
        visible = flag;
    }
    //Do NOT modify this function
    public boolean getVisible() {
        return visible;
    }
    //Do NOT modify this function
    public void addTFChangeListener(TFChangeListener l) {
        if (!listeners.contains(l)) {
            listeners.add(l);
        }
    }
    //Do NOT modify this function
    public abstract void visualize(GL2 gl);
}
