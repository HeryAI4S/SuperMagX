import com.comsol.model.Model;
import com.comsol.model.util.ModelUtil;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Locale;

/** Independent COMSOL 6.4 transient benchmark for a rectangular CuNb ring. */
public class CoupledLossBenchmark {
  static final String DEFAULT_SPEC = "single_a,6e-9,0.02,0.02,0.01,0.001,0.06,0.00025";
  static String units(double v, String unit) {
    return String.format(Locale.US, "%.15g[%s]", v, unit);
  }

  public static void main(String[] args) throws Exception {
    // name rho width height ramp hmax end dt [flatTime]
    // For two independent supplies add flatTime secondR secondZ scale delay.
    if (args == null || args.length == 0)
      args = DEFAULT_SPEC.split(",");
    if (args.length != 8 && args.length != 9 && args.length != 10 && args.length != 13)
      throw new IllegalArgumentException("name rho width height ramp hmax end dt [flat [airExtent or r2 z2 scale delay]]");
    boolean second = args.length == 13;
    String name = args[0];
    double rho = Double.parseDouble(args[1]);
    double width = Double.parseDouble(args[2]);
    double height = Double.parseDouble(args[3]);
    double ramp = Double.parseDouble(args[4]);
    double hmax = Double.parseDouble(args[5]);
    double end = Double.parseDouble(args[6]);
    double dt = Double.parseDouble(args[7]);
    Path out = Paths.get("validation", "comsol", "output");

    Model model = ModelUtil.create("benchmark");
    model.label(name + ".mph");
    model.param().set("rho", units(rho, "ohm*m"));
    model.param().set("w", units(width, "m"));
    model.param().set("h_cond", units(height, "m"));
    model.param().set("Tr", units(ramp, "s"));
    model.param().set("Tp", units(args.length == 9 ? Double.parseDouble(args[8]) : 0.001, "s"));
    model.param().set("Ipk", "1000[A]");
    model.param().set("airR", units(args.length == 10 ?
        Double.parseDouble(args[9]) : 1.0, "m"));
    model.param().set("Iwave", "Ipk*(min(max(t/Tr,0),1)-min(max((t-Tr-Tp)/Tr,0),1))");
    if (second) {
      model.param().set("r2", units(Double.parseDouble(args[9]), "m"));
      model.param().set("z2", units(Double.parseDouble(args[10]), "m"));
      model.param().set("scale2", args[11]);
      model.param().set("delay2", units(Double.parseDouble(args[12]), "s"));
      model.param().set("Iwave2", "Ipk*scale2*(min(max((t-delay2)/Tr,0),1)-min(max((t-delay2-Tr-Tp)/Tr,0),1))");
    }

    model.component().create("comp1", true);
    model.component("comp1").geom().create("geom1", 2);
    model.component("comp1").geom("geom1").axisymmetric(true);
    model.component("comp1").geom("geom1").lengthUnit("m");
    model.component("comp1").geom("geom1").create("air", "Rectangle");
    model.component("comp1").geom("geom1").feature("air")
        .set("pos", new String[]{"0", "-airR"});
    model.component("comp1").geom("geom1").feature("air")
        .set("size", new String[]{"airR", "2*airR"});
    model.component("comp1").geom("geom1").create("metal", "Rectangle");
    model.component("comp1").geom("geom1").feature("metal")
        .set("pos", new String[]{"0.1", "-h_cond/2"});
    model.component("comp1").geom("geom1").feature("metal")
        .set("size", new String[]{"w", "h_cond"});
    model.component("comp1").geom("geom1").feature("metal").set("selresult", true);
    if (second) {
      model.component("comp1").geom("geom1").create("metal2", "Rectangle");
      model.component("comp1").geom("geom1").feature("metal2")
          .set("pos", new String[]{"r2", "z2-h_cond/2"});
      model.component("comp1").geom("geom1").feature("metal2")
          .set("size", new String[]{"w", "h_cond"});
      model.component("comp1").geom("geom1").feature("metal2")
          .set("selresult", true);
    }
    model.component("comp1").geom("geom1").run();
    String metalSel = "geom1_metal_dom";
    if (second) {
      model.component("comp1").selection().create("allMetal", "Union");
      model.component("comp1").selection("allMetal").set("input",
          new String[]{"geom1_metal_dom", "geom1_metal2_dom"});
    }
    String allMetal = second ? "allMetal" : metalSel;
    System.out.println("METAL_DOMAINS=" +
        model.component("comp1").selection(allMetal).entities(2).length);

    model.component("comp1").material().create("matAir", "Common");
    model.component("comp1").material("matAir").label("Vacuum");
    model.component("comp1").material("matAir").propertyGroup("def")
        .set("relpermeability", "1");
    model.component("comp1").material("matAir").propertyGroup("def")
        .set("electricconductivity", "0[S/m]");
    model.component("comp1").material().create("matMetal", "Common");
    model.component("comp1").material("matMetal").selection().named(allMetal);
    model.component("comp1").material("matMetal").label("CuNb constant rho");
    model.component("comp1").material("matMetal").propertyGroup("def")
        .set("relpermeability", "1");
    model.component("comp1").material("matMetal").propertyGroup("def")
        .set("electricconductivity", "1/rho");

    model.component("comp1").physics().create("mf", "InductionCurrents", "geom1");
    model.component("comp1").physics("mf").create("coil1", "Coil", 2);
    model.component("comp1").physics("mf").feature("coil1").selection()
        .named(metalSel);
    model.component("comp1").physics("mf").feature("coil1").set("ICoil", "Iwave");
    if (second) {
      model.component("comp1").physics("mf").create("coil2", "Coil", 2);
      model.component("comp1").physics("mf").feature("coil2").selection()
          .named("geom1_metal2_dom");
      model.component("comp1").physics("mf").feature("coil2")
          .set("ICoil", "Iwave2");
    }

    model.component("comp1").mesh().create("mesh1");
    model.component("comp1").mesh("mesh1").create("ftri1", "FreeTri");
    model.component("comp1").mesh("mesh1").feature("ftri1")
        .create("size1", "Size");
    model.component("comp1").mesh("mesh1").feature("ftri1")
        .feature("size1").selection().named(allMetal);
    model.component("comp1").mesh("mesh1").feature("ftri1")
        .feature("size1").set("custom", true);
    model.component("comp1").mesh("mesh1").feature("ftri1")
        .feature("size1").set("hmax", units(hmax, "m"));
    model.component("comp1").mesh("mesh1").run();
    System.out.println("MESH_BUILT=" + name);

    model.study().create("std1");
    model.study("std1").create("time", "Transient");
    model.study("std1").feature("time").set("tlist",
        String.format(Locale.US, "range(0,%.15g,%.15g)", dt, end));
    model.study("std1").setGenPlots(false);
    model.study("std1").createAutoSequences("all");
    model.sol("sol1").runAll();
    System.out.println("SOLVED=" + name);

    model.result().numerical().create("power", "IntSurface");
    model.result().numerical("power").selection().named(allMetal);
    // Pin the axisymmetric volume convention; never rely on the node default.
    model.result().numerical("power").set("intvolume", "on");
    model.result().numerical("power").set("expr", "mf.Qh");
    model.result().numerical("power").set("unit", "W");
    model.result().numerical("power").set("innerinput", "all");
    double[][] powers = model.result().numerical("power").getReal();
    System.out.println("POWER_SHAPE=" + powers.length + "x" + powers[0].length);
    model.result().numerical().create("powerSection", "IntSurface");
    model.result().numerical("powerSection").selection().named(allMetal);
    model.result().numerical("powerSection").set("intvolume", "off");
    model.result().numerical("powerSection").set("expr", "2*pi*r*mf.Qh");
    model.result().numerical("powerSection").set("unit", "W");
    model.result().numerical("powerSection").set("innerinput", "all");
    double[] sectionPowers = model.result().numerical("powerSection").getReal()[0];
    double axisymError = 0;
    for (int j=0; j<powers[0].length; j++)
      axisymError = Math.max(axisymError,
          Math.abs(powers[0][j]-sectionPowers[j]) /
          Math.max(1.0, Math.abs(powers[0][j])));
    System.out.println(String.format(Locale.US,
        "POWER_AXISYM_MAX_REL=%.15g", axisymError));
    if (axisymError > 1e-9)
      throw new IllegalStateException("axisymmetric integration mismatch");
    model.result().numerical().create("current", "EvalGlobal");
    model.result().numerical("current").set("expr", "mf.ICoil_1");
    model.result().numerical("current").set("innerinput", "all");
    double[][] coilCurrents = model.result().numerical("current").getReal();
    System.out.println("CURRENT_SHAPE=" + coilCurrents.length + "x" + coilCurrents[0].length);
    if (second) {
      model.result().numerical().create("current2", "EvalGlobal");
      model.result().numerical("current2").set("expr", "mf.ICoil_2");
      model.result().numerical("current2").set("innerinput", "all");
      double[][] values2 = model.result().numerical("current2").getReal();
      for (int j=0; j<values2[0].length; j++)
        System.out.println(String.format(Locale.US, "CURRENT2,%.15g,%.15g",
            j*dt, values2[0][j]));
    }
    for (int j=0; j<powers[0].length; j++)
      System.out.println(String.format(Locale.US, "DATA,%.15g,%.15g",
          j*dt, powers[0][j]));
    for (int j=0; j<coilCurrents[0].length; j++)
      System.out.println(String.format(Locale.US, "CURRENT,%.15g,%.15g",
          j*dt, coilCurrents[0][j]));
    String[] diagnostics = {"1", "mf.Jphi/(2*pi*r)",
                            "rho*mf.Jphi^2", "mf.Qh"};
    for (int k=0; k<diagnostics.length; k++) {
      String tag = "diag" + k;
      model.result().numerical().create(tag, "IntSurface");
      model.result().numerical(tag).selection().named(allMetal);
      model.result().numerical(tag).set("intvolume", "on");
      model.result().numerical(tag).set("expr", diagnostics[k]);
      model.result().numerical(tag).set("innerinput", "all");
      double[][] value = model.result().numerical(tag).getReal();
      System.out.println(String.format(Locale.US, "DIAG,%s,%.15g",
          diagnostics[k], value[0][value[0].length-1]));
    }
    model.save(out.resolve(name + ".mph").toString());
    System.out.println("OUTPUT=" + name);
  }
}
