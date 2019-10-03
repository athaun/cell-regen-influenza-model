import java.io.*;

class runTest {

    static double[] runParameters = {0.003, 0.01, 0.03, 0.1, 0.3, 1, 3}; // Mean of exponetial distro for regeneration

    public static void main(String[] args) {




      for (int runCount = 0; runCount < runParameters.length; runCount ++) {
        System.out.println("Running viral transmission Code with regeneration parameter \"" + runParameters[runCount] + "\".");
        runBashCmd("python3 2DProbabilityViralTransmission.py " + runParameters[runCount]); // Run each test individually with one of the parameters
      }

    }

    public static String runBashCmd(String command) {
      /**
        Runs a command in the BASH, then returns the output
      */
      String output = null;
      try {
        // System.out.println("ran command: " + command); // Prints out the command being run
        Process process = Runtime.getRuntime().exec(command); // Run the command in the shell
        BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream())); // Gets output
        while ((output = reader.readLine()) != null) {
          System.out.println(output); // Prints the commands output
        }
      } catch (Exception e) {
        System.out.println("Error running command | Exception caught: " + e);
      }
      return output;
    }
}
