
# Title: compare_runtimes.py

# Description:
# Compile sms-mi-reg registration runtimes from the generated logfile.

# Created on: May 2024
# Created by: Joshua Auger (joshua.auger@childrens.harvard.edu), Computational Radiology Lab, Boston Children's Hospital


import re
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os


def extract_values_from_lines(log_file_path, search_string):
    values = []

    # Create a regex pattern to extract the numerical value following the search string
    pattern = re.compile(rf"{re.escape(search_string)}\s*:\s*([\d\.]+)")

    with open(log_file_path, 'r') as file:
        for line in file:
            match = pattern.search(line)
            if match:
                value = match.group(1)
                # Convert to int if it's a whole number, otherwise to float
                values.append(int(value) if value.isdigit() else float(value))

    return values


def plot_runtimes(df, output_plot_path=None):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

    # Plot runtime for sms-mi-reg and nlopt
    ax1.scatter(df.index, df['sms-mi-reg call runtime (sec)'], color='black', marker='o', s=14, label='sms-mi-reg subprocess call (sec)')
    smsmireg_mean = df['sms-mi-reg call runtime (sec)'].mean()
    ax1.axhline(y=smsmireg_mean, color='black', linestyle='--', label=f'sms-mi-reg subprocess call mean = {smsmireg_mean:.3f} s')
    ax1.scatter(df.index, df['nlopt runtime (sec)'], edgecolor='red', facecolor='none', marker='o', label='nlopt optimizer runtime (sec)')
    nlopt_mean = df['nlopt runtime (sec)'].mean()
    ax1.axhline(y=nlopt_mean, color='red', linestyle='--', label=f'nlopt optimizer mean = {nlopt_mean:.3f} s')
    ax1.set_ylim(bottom=0)
    ax1.set_ylabel('Time (sec)')
    ax1.set_title('Registration Runtimes')
    ax1.legend(loc='upper left')
    ax1.grid(True)

    # Plot cumulative elapsed runtime
    final_cumulative_time = df['cumulative runtime (sec)'].iloc[-1]
    ax2.plot(df.index, df['cumulative runtime (sec)'], color='blue', marker='o', linestyle='-', linewidth=1, markersize=5, label=f'Cumulative runtime: {final_cumulative_time:.3f} sec')
    # Add the reference 1-1 line in red
    ax2.plot(df.index, df.index, color='red', linestyle='-', linewidth=1, label='1 reg/sec reference')
    ax2.set_ylim(bottom=0)
    ax2.set_xlabel('Registration instance')
    ax2.set_ylabel('Elapsed time (sec)')
    ax2.legend(loc='upper left')
    ax2.grid(True)

    # Save the plot if output_plot_path is specified
    if output_plot_path:
        plt.savefig(output_plot_path)
        print(f"Plot saved to {output_plot_path}")

    plt.tight_layout()
    plt.show()


def report_statistics(column):
    mean_val = column.mean()
    stddev_val = column.std()
    min_val = column.min()
    max_val = column.max()
    median_val = column.median()

    print(f"{column.name} stats:")
    print(f" Mean : {mean_val:.4f}")
    print(f" Std Dev : {stddev_val:.4f}")
    print(f" Median : {median_val:.4f}")
    print(f" Min : {min_val:.4f}")
    print(f" Max : {max_val:.4f}")
    print()


def main():
    parser = argparse.ArgumentParser(description='Process a log file to extract and compare runtimes.')
    parser.add_argument('logfile', type=str, help='The path to the log file.')
    args = parser.parse_args()

    log_file_path = args.logfile
    log_file_name = os.path.basename(log_file_path)

    try:
        nlopt_times = extract_values_from_lines(log_file_path, "NLopt elapsed time (sec)")
        nfunc_evals = extract_values_from_lines(log_file_path, "Number of function evals is")
        smsmireg_times = extract_values_from_lines(log_file_path, "Registration call elapsed runtime (sec)")
        total_elapsed_time = extract_values_from_lines(log_file_path, "Total elapsed time (sec)")

        # Ensure lists have the same length
        min_length = min(len(nlopt_times), len(nfunc_evals), len(smsmireg_times), len(total_elapsed_time))
        if min_length == 0:
            print("Error: No valid runtime data found in the log file.")
            return

        # Create DataFrame
        data = {
            'sms-mi-reg call runtime (sec)': smsmireg_times[:min_length],
            'nlopt runtime (sec)': nlopt_times[:min_length],
            'N func evals (n)': nfunc_evals[:min_length],
            'cumulative runtime (sec)': total_elapsed_time[:min_length]
        }
        df = pd.DataFrame(data)

        # Save CSV file
        csv_file_name = f"runtimes_{os.path.splitext(log_file_name)[0]}.csv"
        output_csv_path = os.path.join(os.path.dirname(log_file_path), csv_file_name)
        df.to_csv(output_csv_path, index=False)
        print(f"CSV file saved: {output_csv_path}")

        # Report statistics
        report_statistics(df['sms-mi-reg call runtime (sec)'])
        report_statistics(df['nlopt runtime (sec)'])
        report_statistics(df['N func evals (n)'])

        # Generate plot
        plot_file_name = f"runtimes_{os.path.splitext(log_file_name)[0]}.png"
        output_plot_path = os.path.join(os.path.dirname(log_file_path), plot_file_name)
        plot_runtimes(df, output_plot_path)

    except FileNotFoundError:
        print(f"Error: Log file '{log_file_path}' not found.")
    except pd.errors.EmptyDataError:
        print(f"Error: Extracted data resulted in an empty DataFrame.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


if __name__ == "__main__":
    main()
