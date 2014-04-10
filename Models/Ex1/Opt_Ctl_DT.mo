within ;
package Opt_Ctl_DT
  "Models for optimal control using differential transformations method"
  model Ex_1 "Example 1: Source of this example in the info layer"

  // State start values
      parameter Real x_0 = 1;
  // The states
     Real x(start = x_0);

     // The control signal
    Modelica.Blocks.Interfaces.RealInput u
      annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));
  equation
     der(x) = -x + u;

    annotation (
      Documentation(info="<HTML>
<p>
The model is a modelica implementation of example-3 from:</p>
<p>
Hwang I, Li J, Du D. Differential Transformation and Its Application to Nonlinear Optimal Control. J. Dyn. Sys., Meas., Control. 2009;131(5):051010-051010-11. 
</p>
<p>
This model is set-up for an infinite time optimal control problem whose optimization class is set up in:./OPT.mop in the class Ex1_Opt
</p>
</HTML>"),
        experiment(StopTime=10),
      Icon(graphics={Line(
            points={{40,4},{40,0},{38,-4},{34,-10},{28,-16},{22,-20},{16,-24},{0,-34},
                {-6,-36},{-10,-36},{-14,-34},{-18,-32},{-22,-28},{-26,-22},{-30,-18},
                {-32,-14},{-36,-8},{-36,-4},{-36,2},{-36,8},{-34,16},{-30,20},{-24,
                26},{-14,32},{2,42},{12,50},{20,54},{30,56},{40,56},{46,54},{54,48},
                {60,42},{64,34},{66,28},{66,22},{66,18},{66,12},{64,6},{62,2},{58,
                -2},{48,-12},{38,-20},{24,-30},{6,-40},{-4,-44},{-14,-44},{-22,-42},
                {-30,-36},{-36,-30},{-42,-24},{-46,-18},{-50,-10},{-50,0},{-48,12},
                {-44,20},{-36,26},{-30,30}},
            color={0,0,255},
            smooth=Smooth.None)}));
  end Ex_1;
  annotation (uses(Modelica(version="3.2")));
end Opt_Ctl_DT;
