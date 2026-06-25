import sys

def calculate_average(numbers):
    return sum(numbers) / len(numbers)

if __name__ == "__main__":
    # Read numbers from stdin
    input_data = sys.stdin.read()

    # Convert input data into a list of floats
    numbers = [float(num) for num in input_data.split()]

    # Calculate the average
    average = calculate_average(numbers)

    # Print the result
    print(f"The average of the array is: {average}")
