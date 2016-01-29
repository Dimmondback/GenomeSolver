
import java.awt.Color;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.ArrayList;
import java.util.List;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.JTextPane;
import javax.swing.border.Border;
import javax.swing.border.LineBorder;

/**
 *
 * @author Nick Saric // Date: Oct 01, 2015 // Used for Bioinformatics
 * similarity table parsing and Genome comparisons. Can be done side by side or
 * substring Genome comparison.
 */
public class GenomicSimilaritySolver {

    public static ArrayList<String> tableLetters;
    public static List<Integer> inputTableArray;
    static JFrame frame1 = new JFrame("Genomic Similarity Solver");
    static JFrame frame2 = new JFrame("Needleman-Wunsch Solver");
    static int whichFrame = 1;
    static String last1 = "";
    static String last2 = "";
    static int penalty = 0; //How many gaps were used. Gaps can be of any length.
    //calculated_penalty = -init_gapscore - (#_moves - 1){E::(-1)} #gaps starts at 1. [As gap gets larger, penalty gets less!]

    public static int posGet(int row, int column, List<Integer> array, int length) {
        int pos = (column * length) + row;
        int retThis = array.get(pos);
        return retThis;
    }

    public static int[][] matrixGen(String input1, String input2, List<Integer> table, String gapScore, boolean smith) {
        int gapScoreInt = Integer.parseInt(gapScore);
        int[][] retThis = new int[input1.length() + 1][input2.length() + 1];
        for (int i = 0; i < input1.length() + 1; i++) {
            for (int j = 0; j < input2.length() + 1; j++) {
                if (i == 0 || j == 0) {
                    if (j == 0) {
                        if (smith) {
                            for (int k = 0; k < i; k++) {
                                int temp = Integer.parseInt(gapScore) - (k - 1) * (-1);
                                if (true) {
                                    retThis[i][j] = temp;
                                }
                            }
                        } else {
                            retThis[i][j] = gapScoreInt * j;
                        }
                    } else if (i == 0) {
                        if (smith) {
                            for (int k = 0; k < j; k++) {
                                int temp = Integer.parseInt(gapScore) - (k - 1) * (-1);
                                if (true) {
                                    retThis[i][j] = temp;
                                }
                            }
                        } else {
                            retThis[i][j] = gapScoreInt * j;
                        }
                        retThis[i][j] = gapScoreInt * i;
                    } else {
                        //Comment
                    }
                }
            }
        }
        return retThis;
    }

    public static int matrixFill() {
        //echh
        return 1;
    }

    public static int smitWat(ArrayList<String> letterArrayList, List<Integer> inputNums, String[] input1, String[] input2) {
        int x = 0;
        int y = 0;
        int retThis = 0;

        if (input1.length <= input2.length) {
            for (int i = 0; i < input1.length; i++) {
                x = letterArrayList.indexOf(input1[i]);
                last1 = input1[i];
                y = letterArrayList.indexOf(input2[i]);
                last2 = input2[i];
                if (last1.equals("") || last2.equals("")) {
                    continue;
                }
                retThis = retThis + posGet(x, y, inputNums, letterArrayList.size());
            }
        } else {
            retThis = sideBySide(letterArrayList, inputNums, input2, input1);
        }
        return retThis;
    }

    public static int sideBySide(ArrayList<String> letterArrayList, List<Integer> inputNums, String[] input1, String[] input2) {
        int x = 0;
        int y = 0;
        int retThis = 0;

        if (input1.length <= input2.length) {
            for (int i = 0; i < input1.length; i++) {
                x = letterArrayList.indexOf(input1[i]);
                last1 = input1[i];
                y = letterArrayList.indexOf(input2[i]);
                last2 = input2[i];
                if (last1.equals("") || last2.equals("")) {
                    continue;
                }
                retThis = retThis + posGet(x, y, inputNums, letterArrayList.size());
            }
        } else {
            retThis = sideBySide(letterArrayList, inputNums, input2, input1);
        }
        return retThis;
    }

    public static int needWun(ArrayList<String> letterArrayList, List<Integer> inputNums, String[] input1, String[] input2) {
        int x = 0;
        int y = 0;
        int retThis = 0;

        if (input1.length <= input2.length) {
            for (int i = 0; i < input1.length; i++) {
                x = letterArrayList.indexOf(input1[i]);
                last1 = input1[i];
                y = letterArrayList.indexOf(input2[i]);
                last2 = input2[i];
                if (last1.equals("") || last2.equals("")) {
                    continue;
                }
                retThis = retThis + posGet(x, y, inputNums, letterArrayList.size());
            }
        } else {
            retThis = sideBySide(letterArrayList, inputNums, input2, input1);
        }
        return retThis;
    }

    public static int subInt(ArrayList<String> letterArrayList, List<Integer> inputNums, String[] input1, String[] input2) {
        int retThis = 0;

        if (input1.length <= input2.length) {
            for (int i = 0; i < input1.length; i++) {
                for (int j = 0; j < input2.length; j++) {
                    try {
                        if (input1[i].equals("") || input2[j].equals("")) {
                            if (input1[i].equals("")) {
                                break;
                            } else {
                                continue;
                            }
                        }
                        if (input1[i].equals(input2[j])) {
                            ArrayList<String> sub1 = new ArrayList<String>();
                            for (int k = 0; k < input1.length; k++) {
                                if (!input1[k].equals("")) {
                                    sub1.add(input1[k]);
                                }
                            }
                            ArrayList<String> sub2 = new ArrayList<String>();
                            for (int k = 0; k < sub1.size(); k++) {
                                if (input2[0].equals("")) {
                                    sub2.add(input2[j - i + k + 1]);
                                } else {
                                    sub2.add(input2[j - i + k]);
                                }
                            }
                            int x = letterArrayList.indexOf(input1[i]);
                            last1 = input1[i];
                            int y = letterArrayList.indexOf(input2[j]);
                            last2 = input2[j];
                            int test = posGet(x, y, inputNums, letterArrayList.size());
                            if (test > retThis) {
                                retThis = test;
                            }
                        }
                    } catch (IndexOutOfBoundsException e) {
                    }
                }
            }
        } else {
            retThis = subInt(letterArrayList, inputNums, input2, input1);
        }
        return retThis;
    }

    public static List<Integer> tableParse(String input) {
        input = input.replaceAll("\\s", ",");
        String[] split = input.split(",");
        List<Integer> nums = new ArrayList<Integer>();
        ArrayList<String> letters = new ArrayList<String>();
        for (int i = 0; i < split.length; i++) {
            if (java.util.regex.Pattern.matches("-?\\d+", split[i])) {
                nums.add(Integer.parseInt(split[i].replaceAll("\\s", "")));
            } else if (split[i].replaceAll("\\s", "").length() > 0) {
                if (!letters.contains(split[i].replaceAll("\\s", ""))) {
                    letters.add(split[i].replaceAll("\\s", ""));
                }
            }
        }
        tableLetters = letters;
        return nums;
    }

    public static void windowSetup() {
        final JTextArea arrayText = new JTextArea();
        final JTextField textArea1 = new JTextField();
        final JTextField textArea2 = new JTextField();
        final JTextPane resultArea = new JTextPane();
        final JCheckBox checkbox = new JCheckBox("Check this box for substring comparison method.");
        JButton submit = new JButton("Submit");
        Border border = new LineBorder(Color.GRAY);
        frame1.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
        frame1.addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                exitWindow();
            }
        });
        frame1.setResizable(false);
        frame1.setLayout(null);
        frame1.setBounds((int) Toolkit.getDefaultToolkit().getScreenSize().getWidth() / 8 - frame1.getWidth() / 2, (int) Toolkit.getDefaultToolkit().getScreenSize().getHeight() / 8 - frame1.getHeight() / 2, 900, 500);
        arrayText.setBounds(1, 1, frame1.getWidth() / 2 - 3, frame1.getHeight() - 91);
        arrayText.setBorder(border);
        arrayText.setLineWrap(true);
        arrayText.setWrapStyleWord(false);
        resultArea.setBounds(arrayText.getWidth() + 2, 1, frame1.getWidth() / 2 - 7, arrayText.getHeight() - 31);
        resultArea.setBorder(border);
        resultArea.setFocusable(false);
        resultArea.setText("Please paste your table of elements to the left.\n\n\nPlease paste your two genomes to compare below.\n\n\nYou can switch to the Needleman-Wunsch Solver by hitting the [X] in the top left of the window and clicking on the \"Switch\" button");
        submit.setBounds(resultArea.getX() + resultArea.getWidth() - resultArea.getWidth() / 4, resultArea.getHeight() + 2, resultArea.getWidth() / 4, 30);
        submit.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                inputTableArray = tableParse(arrayText.getText());
                String rawIn1 = textArea1.getText().replaceAll("\\s", "");
                String rawIn2 = textArea2.getText().replaceAll("\\s", "");
                String[] input1 = textArea1.getText().replaceAll("\\s", "").split("");
                String[] input2 = textArea2.getText().replaceAll("\\s", "").split("");
                if (inputTableArray.size() == 0 || rawIn1.length() == 0 || rawIn2.length() == 0) {
                    resultArea.setText("One or more of your inputs were empty.\n\n\nPlease ensure all fields are filled.");
                    return;
                }
                try {
                    if (checkbox.isSelected()) {
                        //Substring method.
                        resultArea.setText("Sub final similarity is: " + subInt(tableLetters, inputTableArray, input1, input2));
                    } else {
                        //Side by side comparison method.
                        resultArea.setText("Side by side final similarity is: " + sideBySide(tableLetters, inputTableArray, input1, input2));
                    }
                } catch (IndexOutOfBoundsException z) {
                    resultArea.setText("Your table cannot compare the two letters " + last1 + " and " + last2 + ".\n\n\nPlease ensure your table is square and completely filled.");
                }
            }
        });
        checkbox.setBounds(resultArea.getX(), resultArea.getHeight() + 1, resultArea.getWidth() - submit.getWidth() - 1, 30);
        textArea1.setBounds(1, arrayText.getHeight() + 2, frame1.getWidth() - 7, 30);
        textArea2.setBounds(1, textArea1.getY() + textArea1.getHeight() + 0, frame1.getWidth() - 7, 30);
        frame1.add(arrayText);
        frame1.add(textArea1);
        frame1.add(textArea2);
        frame1.add(resultArea);
        frame1.add(checkbox);
        frame1.add(submit);
        frame1.setVisible(true);
    }

    public static void windowSetup2() {
        final JTextArea arrayText = new JTextArea();
        final JTextField textArea1 = new JTextField();
        final JTextField textArea2 = new JTextField();
        final JTextPane resultArea = new JTextPane();
        final JTextField gapScore = new JTextField();
        final JTextField extScore = new JTextField();
        final JTextField gapScoreText = new JTextField();
        final JTextField extScoreText = new JTextField();
        final JCheckBox checkbox = new JCheckBox("Smith-Waterman version.");
        JButton submit = new JButton("Submit");
        Border border = new LineBorder(Color.GRAY);
        frame2.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
        frame2.addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                exitWindow();
            }
        });
        frame2.setResizable(false);
        frame2.setLayout(null);
        frame2.setBounds((int) Toolkit.getDefaultToolkit().getScreenSize().getWidth() / 8 - frame2.getWidth() / 2, (int) Toolkit.getDefaultToolkit().getScreenSize().getHeight() / 8 - frame2.getHeight() / 2, 900, 500);
        arrayText.setBounds(1, 1, frame2.getWidth() / 2 - 3, frame2.getHeight() - 91);
        arrayText.setBorder(border);
        arrayText.setLineWrap(true);
        arrayText.setWrapStyleWord(false);
        resultArea.setBounds(arrayText.getWidth() + 2, 1, frame2.getWidth() / 2 - 7, arrayText.getHeight() - 31 * 2);
        resultArea.setBorder(border);
        resultArea.setFocusable(false);
        resultArea.setText("Please paste your table of elements to the left.\n\n\nPlease paste your two genomes to compare below.\n\n\nYou can switch to the Genomic Similarity Solver by hitting the [X] in the top left of the window and clicking on the \"Switch\" button");
        submit.setBounds(resultArea.getX() + resultArea.getWidth() - resultArea.getWidth() / 4, resultArea.getHeight() + 2, resultArea.getWidth() / 4, 60);
        //Fix methods!
        submit.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                inputTableArray = tableParse(arrayText.getText());
                String rawIn1 = textArea1.getText().replaceAll("\\s", "");
                String rawIn2 = textArea2.getText().replaceAll("\\s", "");
                String[] input1 = textArea1.getText().replaceAll("\\s", "").split("");
                String[] input2 = textArea2.getText().replaceAll("\\s", "").split("");
                if (inputTableArray.size() == 0 || rawIn1.length() == 0 || rawIn2.length() == 0 || gapScore.getText().length() == 0) {
                    resultArea.setText("One or more of your inputs were empty.\n\n\nPlease ensure all fields are filled.");
                    return;
                }
                try {
                    if (checkbox.isSelected()) {
                        //Smith-Waterman method.
                        resultArea.setText("Smith-Waterman score is: " + smitWat(tableLetters, inputTableArray, input1, input2));
                    } else {
                        //Needleman-Wunsch method.
                        resultArea.setText("Needleman-Wunsch final score is: " + needWun(tableLetters, inputTableArray, input1, input2));
                    }
                } catch (IndexOutOfBoundsException z) {
                    resultArea.setText("Your table cannot compare the two letters " + last1 + " and " + last2 + ".\n\n\nPlease ensure your table is square and completely filled.");
                }
            }
        });
        gapScore.setBounds(resultArea.getX(), resultArea.getHeight() + 2, 30, 31);
        extScore.setBounds(resultArea.getX(), gapScore.getY() + gapScore.getHeight() + 2, 30, 31);
        gapScoreText.setBounds(gapScore.getX() + gapScore.getWidth() + 1, resultArea.getHeight() + 2, 125, 31);
        gapScoreText.setText("Gap score goes here.");
        extScoreText.setBounds(gapScore.getX() + gapScore.getWidth() + 1, resultArea.getHeight() + 2, 125, 31);
        extScoreText.setText("Extension penalty goes here.");
        checkbox.setLocation(gapScoreText.getX() + gapScoreText.getWidth() + 1, gapScoreText.getY() + gapScoreText.getHeight() + 1);
        checkbox.setSize(submit.getX() - checkbox.getX() - 1, 30);
        gapScoreText.setEditable(false);
        textArea1.setBounds(1, arrayText.getHeight() + 2, frame2.getWidth() - 7, 30);
        textArea2.setBounds(1, textArea1.getY() + textArea1.getHeight() + 0, frame2.getWidth() - 7, 30);
        frame2.add(arrayText);
        frame2.add(textArea1);
        frame2.add(textArea2);
        frame2.add(resultArea);
        frame2.add(gapScore);
        frame2.add(submit);
        frame2.add(checkbox);
        frame2.add(gapScoreText);
        frame2.setVisible(false);
    }

    public static void exitWindow() {
        final JFrame exitWindow = new JFrame("Exit");
        JButton switcher = new JButton("Switch");
        JButton yes = new JButton("Exit");
        JButton no = new JButton("Cancel");
        switcher.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                if (whichFrame == 1) {
                    whichFrame = 2;
                    frame1.setVisible(false);
                    frame1.setFocusable(false);
                    frame2.setLocation(frame1.getLocation());
                    frame2.setFocusable(true);
                    frame2.setVisible(true);
                    frame2.setFocusableWindowState(true);
                    exitWindow.dispose();
                } else {
                    whichFrame = 1;
                    frame2.setVisible(false);
                    frame2.setFocusable(false);
                    frame1.setLocation(frame2.getLocation());
                    frame1.setFocusable(true);
                    frame1.setVisible(true);
                    frame1.setFocusableWindowState(true);
                    exitWindow.dispose();
                }
            }
        });
        yes.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                System.exit(0);
            }
        });
        no.addActionListener(new ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                if (whichFrame == 1) {
                    frame1.setFocusableWindowState(true);
                } else {
                    frame2.setFocusableWindowState(true);
                }
                exitWindow.dispose();
            }
        });
        switcher.setSize(100, 40);
        yes.setSize(100, 40);
        no.setSize(100, 40);
        exitWindow.setSize(380, 85);
        exitWindow.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        exitWindow.setUndecorated(true);
        exitWindow.getRootPane().setBorder(new LineBorder(Color.getHSBColor((float) .6, (float) .1, (float) .9), 3));
        if (whichFrame == 1) {
            exitWindow.setLocation(frame1.getX() + frame1.getWidth() / 2 - exitWindow.getWidth() / 2, frame1.getY() + frame1.getHeight() / 2 - exitWindow.getHeight() / 2);
        } else {
            exitWindow.setLocation(frame2.getX() + frame2.getWidth() / 2 - exitWindow.getWidth() / 2, frame2.getY() + frame2.getHeight() / 2 - exitWindow.getHeight() / 2);
        }
        if (whichFrame == 1) {
            frame1.setFocusableWindowState(false);
        } else {
            frame2.setFocusableWindowState(false);
        }
        exitWindow.setLayout(null);
        exitWindow.setResizable(false);
        exitWindow.add(switcher);
        exitWindow.add(yes);
        exitWindow.add(no);
        if (whichFrame == 1) {
            frame1.setFocusable(false);
        } else {
            frame2.setFocusable(false);
        }
        switcher.setLocation(20, 20);
        yes.setLocation(switcher.getX() + switcher.getWidth() + 20, 20);
        no.setLocation(yes.getX() + yes.getWidth() + 20, 20);
        exitWindow.setVisible(true);
    }

    public static void main(String[] args) {
        windowSetup();
        windowSetup2();
    }
}