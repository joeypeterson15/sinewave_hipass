#include <iostream>
#include <cmath>
#include <fstream>
#include <array>
#include <complex>

using namespace std;

const int sampleRate = 44100;
const int bitDepth = 16;
const int duration = 2;
const int sampleLength = sampleRate * duration;
const float nyquisFrequency = sampleRate / 2;
const float c4NoteHz = 262;
const float a4NoteHz = 440;
const float e5NoteHz = 659.25;
const float g5NoteHz = 783.99;
const double e = 2.71828182845904523536;


class SineOscillator {
    float frequency, amplitude, angle = 0.0f, offset = 0.0f;
public:
    SineOscillator(float freq, float amp) : frequency(freq), amplitude(amp) {
        offset = 2 * M_PI * frequency / sampleRate;
    }
    float process() {
        auto sample = amplitude * sin(angle);
        angle += offset;
        return sample;
        // Asin(2pif/sr)
    }
};

void writeToFile(ofstream &file, int value, int size) {
    file.write(reinterpret_cast<const char*> (&value), size);
}

int main() {
    int N = sampleRate * duration;
    SineOscillator sineOscillatorA(a4NoteHz, 0.125);
    SineOscillator sineOscillatorC(c4NoteHz, 0.125);
    SineOscillator sineOscillatorE(e5NoteHz, 0.125);
    SineOscillator sineOscillatorG(g5NoteHz, 0.125);

    // collect sound data
    std::array <float, sampleLength> samples;
    for (int i = 0; i < N; ++i) {
        float sample = sineOscillatorA.process() + sineOscillatorC.process() + sineOscillatorE.process() + sineOscillatorG.process();
        samples[i] = sample;
    }

    // perform DFT on data:
        std::array<complex<float>, sampleLength> frequencyDomain;
        for (int k = 0; k < N; ++k) {
            complex<float> amplitudeCoefficient = 0;
            for (int n = 0; n < N; ++n) {
                // angle == theta from e^i(theta)
                float angle = -2 * M_PI * k * n / N;
                amplitudeCoefficient += samples[n] * complex<float>(cos(angle), sin(angle));
            }
            // if ( k < 10) {
            //     cout << " First kth amplitudecoeffecient" << amplitudeCoefficient << endl;
            // }
            frequencyDomain[k] = amplitudeCoefficient;
        }

    // process the frequency domain for low pass filter at cutoff frequency = 300 Hz
    // int cutoffBin = cutoffFrequency * sampleRate / N;
    int cutoffFrequency = 300;
    for (int k = 0; k < N / 2; ++k) {
        int freqBin = k * sampleRate / N;
        if (freqBin > cutoffFrequency) {
            frequencyDomain[k] = 0.0;
        }
    }

    for (int k = N / 2 + 1; k < N; ++k) {
        int freqBin = (k - N) * sampleRate / N;  // Negative frequency
        if (freqBin < -cutoffFrequency)  // Now comparing directly with the negative cutoff
            frequencyDomain[k] = 0.0;
    }

    // Perform IDFT
    for (int n = 0; n < N; ++n) {
        complex<float> sum = 0;
        for (int k = 0; k < N; ++k) {
            float angle = 2 * M_PI * k * n / N;
            sum += frequencyDomain[k] * complex<float>(cos(angle), sin(angle));
        }
        samples[n] = sum.real() / float(N);  // Normalize
    }

    ofstream audioFile;
    audioFile.open("waveformAmaj7LoPassCheck.wav", ios::binary);
    //Header chunk
    audioFile << "RIFF";
    audioFile << "----";
    audioFile << "WAVE";

    // Format chunk
    audioFile << "fmt ";
    writeToFile(audioFile, 16, 4); // Size
    writeToFile(audioFile, 1, 2); // Compression code
    writeToFile(audioFile, 1, 2); // Number of channels
    writeToFile(audioFile, sampleRate, 4); // Sample rate
    writeToFile(audioFile, sampleRate * bitDepth / 8, 4 ); // Byte rate
    writeToFile(audioFile, bitDepth / 8, 2); // Block align
    writeToFile(audioFile, bitDepth, 2); // Bit depth

    //Data chunk
    audioFile << "data";
    audioFile << "----";

    int preAudioPosition = audioFile.tellp();

    auto maxAmplitude = pow(2, bitDepth - 1) - 1;
    for(int i = 0; i < N; i++ ) {
        auto sample = samples[i];
        int intSample = static_cast<int> (sample * maxAmplitude);
        writeToFile(audioFile, intSample, 2);
    }
    int postAudioPosition = audioFile.tellp();

    audioFile.seekp(preAudioPosition - 4);
    writeToFile(audioFile, postAudioPosition - preAudioPosition, 4);

    audioFile.seekp(4, ios::beg);
    writeToFile(audioFile, postAudioPosition - 8, 4);

    audioFile.close();
    return 0;
}