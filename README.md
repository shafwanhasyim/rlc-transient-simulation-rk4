# RLC Transient Simulation using Runge-Kutta 4th Order Method (RK4)

## Anggota Kelompok 13
- Aisya Rivelia Azzahra (2306161864)
- Fadhlureza Sebastian (2306161971)
- M. Ikhsan Kurniawan (2306210784)
- Muhammad Raditya Alif Nugroho (2306212745)
- Shafwan Hasyim (2306209113)

## 🧾 Deskripsi Program
Program ini mensimulasikan respons arus dan muatan listrik \( q(t) \) pada rangkaian RLC seri yang diberi eksitasi tegangan sinusoidal. Penyelesaian dilakukan dengan metode numerik **Runge-Kutta orde 4 (RK4)**, dan hasilnya dibandingkan dengan **solusi analitik** yang tersedia dalam buku *Chapman* s(Persamaan 28.11). Perbandingan antara kedua solusi ditampilkan dalam bentuk tabel dan grafik untuk menilai akurasi numerik.

Program ditulis menggunakan bahasa **C** dan menghasilkan output berupa:
- Nilai waktu \( t \)
- Nilai muatan hasil RK4
- Nilai muatan hasil analitiks
- Selisih (error) antara keduanya

