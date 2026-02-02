import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
from math import factorial, ceil

# ==============================
# MODEL M/M/m
# ==============================
def mm_m_full(lam, mu, m, max_p=5):
    lam = float(lam)
    mu = float(mu)
    m = int(m)
    if m < 1:
        raise ValueError("Liczba kanałów musi być >= 1")

    rho = lam / (m * mu)
    if rho >= 1:
        return [
            ("ρ", rho), ("P0", 0), ("Lq", np.inf), ("L", np.inf),
            ("Wq", np.inf), ("W", np.inf), ("m", m),
            ("m_free", 0), ("m_busy", m)
        ] + [("", "")] + [(f"P{k}", np.nan) for k in range(max_p + 1)]

    a = lam / mu

    # P0
    sum_terms = sum((a ** n) / factorial(n) for n in range(m))
    last_term = (a ** m) / (factorial(m) * (1 - rho))
    P0 = 1 / (sum_terms + last_term)

    # Lq wg Erlang C
    Lq = (P0 * (a ** m) * rho) / (factorial(m) * (1 - rho) ** 2)
    Wq = Lq / lam if lam > 0 else 0
    L = Lq + a
    W = Wq + 1 / mu

    # Prawdopodobieństwa P0..Pk
    probs = []
    for k in range(max_p + 1):
        if k < m:
            pk = (a ** k) / factorial(k) * P0
        else:
            pk = (a ** k) / (factorial(m) * m ** (k - m)) * P0
        probs.append((f"P{k}", pk))

    m_busy = a
    m_free = m - m_busy

    metrics = [
        ("ρ", rho),
        ("P0", P0),
        ("Lq", Lq),
        ("L", L),
        ("Wq", Wq),
        ("W", W),
        ("m", float(m)),
        ("m_free", m_free),
        ("m_busy", m_busy),
        ("", ""),
    ] + probs

    return metrics

# ==============================
# OPTYMALIZACJA Z KOSZTAMI C1 i C2
# ==============================
def cockroach_optimization(lam, mu, m_max, c1, c2, n=30, iters=150, step=0.5):
    a = lam / mu
    m_min = int(ceil(a)) if a != int(a) else int(a) + 1
    if m_max < m_min:
        m_max = m_min + 10

    pos = np.random.uniform(m_min, m_max, n)
    best_cost = np.inf
    best_m = m_min

    for _ in range(iters):
        for i in range(n):
            m_try = max(m_min, int(round(pos[i])))
            data = mm_m_full(lam, mu, m_try)
            Wq = [v for k, v in data if k == "Wq"][0]

            # Całkowity koszt: C1 * (średnia liczba klientów w systemie) + C2 * m
            # Średnia liczba klientów w systemie = L = Lq + a
            L = [v for k, v in data if k == "L"][0]
            total_cost = c1 * L + c2 * m_try

            if total_cost < best_cost:
                best_cost = total_cost
                best_m = m_try

        # Aktualizacja pozycji karaluchów
        for i in range(n):
            if np.random.rand() < 0.7:
                pos[i] += step * (best_m - pos[i])
            else:
                pos[i] += step * (2 * np.random.rand() - 1)
            pos[i] = np.clip(pos[i], m_min, m_max)

    return best_m, best_cost

# Mapowanie symboli na pełne nazwy
param_names = {
    "ρ": "Współczynnik obciążenia (ρ)",
    "P0": "Prawdopodobieństwo stanu pustego (P0)",
    "Lq": "Średnia długość kolejki (Lq)",
    "L": "Średnia liczba zgłoszeń w systemie (L)",
    "Wq": "Średni czas oczekiwania w kolejce (Wq)",
    "W": "Średni czas oczekiwania w systemie (W)",
    "m": "Liczba kanałów obsługi (m)",
    "m_free": "Średnia liczba kanałów wolnych (m_free)",
    "m_busy": "Średnia liczba kanałów zajętych (m_busy)",
}

def get_full_name(key):
    if key.startswith("P") and len(key) > 1 and key[1:].isdigit():
        return f"Prawdopodobieństwo {key} zgłoszeń w systemie ({key})"
    return param_names.get(key, key)

# ==============================
# GUI
# ==============================
class QueueApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("System M/M/m/FIFO/∞  - Analiza i Optymalizacja")
        self.geometry("1100x800")

        self.style = ttk.Style(self)
        self.style.theme_use("clam")
        self.style.configure("Treeview", font=("Segoe UI", 10), rowheight=24)
        self.style.configure("Treeview.Heading", font=("Segoe UI", 11, "bold"))
        self.style.configure("TLabel", font=("Segoe UI", 11))

        self.results_before = []
        self.results_after = []

        self.create_header()
        self.create_tabs()

        # Re-wrap long summary texts when resizing window
        self.bind('<Configure>', self._update_wraplength)



    def _update_wraplength(self, event=None):
        """Keep summary labels readable by wrapping text to current window width."""
        try:
            w = max(600, self.winfo_width() - 60)
            # Compare tab summary
            if hasattr(self, "summary_label"):
                self.summary_label.configure(wraplength=w)
            # Optimization tab summary labels (if present)
            if hasattr(self, "objective_label_opt_widget"):
                self.objective_label_opt_widget.configure(wraplength=w)
            if hasattr(self, "efficiency_label_opt_widget"):
                self.efficiency_label_opt_widget.configure(wraplength=w)
        except Exception:
            pass

# --- helpers ---
    def _results_to_dict(self, results):
        """Convert list of (key, value) pairs into a dict, skipping separators."""
        d = {}
        for k, v in results:
            if k == "" and v == "":
                continue
            d[k] = v
        return d

    def _fmt(self, v, digits=4):
        if v is None:
            return "-"
        try:
            if np.isinf(v):
                return "∞"
        except Exception:
            pass
        if isinstance(v, (int, np.integer)):
            return str(int(v))
        if isinstance(v, (float, np.floating)):
            return f"{float(v):.{digits}f}"
        return str(v)

    def create_header(self):
        header = tk.Frame(self, bg="#0f766e", height=60)
        header.pack(fill="x")
        tk.Label(header, text="System Kolejkowy M/M/m/FIFO/∞ \nAnaliza i Optymalizacja z kosztami C1 i C2 (CSO)",
                 bg="#0f766e", fg="white", font=("Segoe UI", 16, "bold")).pack(expand=True)

    def create_tabs(self):
        notebook = ttk.Notebook(self)
        notebook.pack(fill="both", expand=True, padx=10, pady=10)

        self.tab_analysis = ttk.Frame(notebook)
        self.tab_opt = ttk.Frame(notebook)
        self.tab_compare = ttk.Frame(notebook)

        notebook.add(self.tab_analysis, text="Analiza")
        notebook.add(self.tab_opt, text="Optymalizacja")
        notebook.add(self.tab_compare, text="Porównanie")

        self.analysis_tab()
        self.optimization_tab()
        self.compare_tab()

    def analysis_tab(self):
        input_frame = ttk.LabelFrame(self.tab_analysis, text="Parametry systemu", padding=10)
        input_frame.pack(side="left", fill="y", padx=10, pady=10)

        self.lam_var = tk.DoubleVar(value=0.244)
        self.mu_var = tk.DoubleVar(value=0.300)
        self.m_var = tk.IntVar(value=2)
        self.max_p_var = tk.IntVar(value=5)

        # DODANE: Zmienne dla kosztów C1 i C2 na pierwszej stronie
        self.c1_analysis_var = tk.DoubleVar(value=20.0)
        self.c2_analysis_var = tk.DoubleVar(value=50.0)

        for label, var in [("λ (przyjścia/h)", self.lam_var),
                           ("μ (obsługi/h)", self.mu_var),
                           ("m (kanały początkowe)", self.m_var),
                           ("Liczba Pk do wyświetlenia", self.max_p_var),
                           ("Koszt C1 (klient w syst.)", self.c1_analysis_var),  # Nowe pole
                           ("Koszt C2 (kanał)", self.c2_analysis_var)]:  # Nowe pole
            row = tk.Frame(input_frame)
            row.pack(fill="x", pady=4)
            tk.Label(row, text=label, width=28, anchor="w").pack(side="left")
            tk.Entry(row, textvariable=var, width=12).pack(side="left")

        ttk.Button(input_frame, text="Analizuj system", command=self.run_analysis).pack(pady=15)

        # DODANE: Wyświetlanie funkcji celu na dole panelu bocznego
        self.analysis_cost_label = tk.StringVar(value="Koszt J(m): -")
        tk.Label(input_frame, text="Funkcja celu:", font=("Segoe UI", 10, "bold")).pack(pady=(10, 0))
        tk.Label(input_frame, text="J(m) = C1·L + C2·m", font=("Segoe UI", 9, "italic")).pack()
        tk.Label(input_frame, textvariable=self.analysis_cost_label, font=("Segoe UI", 11, "bold"), fg="darkred").pack(
            pady=5)

        self.table_analysis = ttk.Treeview(self.tab_analysis, columns=("param", "val"), show="headings", height=28)
        self.table_analysis.heading("param", text="Parametr")
        self.table_analysis.heading("val", text="Wartość")
        self.table_analysis.column("param", width=450, anchor="w")
        self.table_analysis.column("val", width=180, anchor="center")
        self.table_analysis.pack(side="right", expand=True, fill="both", padx=10, pady=10)

    def run_analysis(self):
        try:
            lam = self.lam_var.get()
            mu = self.mu_var.get()
            m = self.m_var.get()
            max_p = self.max_p_var.get()
            c1 = self.c1_analysis_var.get()
            c2 = self.c2_analysis_var.get()

            self.results_before = mm_m_full(lam, mu, m, max_p)

            self.table_analysis.delete(*self.table_analysis.get_children())
            L_val = 0
            for k, v in self.results_before:
                if k == "L": L_val = v  # Pobieramy L do kosztów

                if k == "" and v == "":
                    self.table_analysis.insert("", "end", values=("", ""), tags=("sep",))
                else:
                    full_param = get_full_name(k)
                    val_str = f"{v:.10f}" if isinstance(v, float) and not np.isinf(v) else "∞"
                    self.table_analysis.insert("", "end", values=(full_param, val_str))

            # OBLICZENIE KOSZTU NA 1 STRONIE
            if not np.isinf(L_val):
                total_cost = c1 * L_val + c2 * m
                self.analysis_cost_label.set(f"Koszt J(m): {total_cost:.4f}")
            else:
                self.analysis_cost_label.set("Koszt J(m): ∞ (niestabilny)")

            self.table_analysis.tag_configure("sep", background="#ebebeb")

        except Exception as e:
            messagebox.showerror("Błąd", str(e))

    def optimization_tab(self):

        # === KONFIGURACJA SIATKI GŁÓWNEJ ===
        self.tab_opt.columnconfigure(0, weight=0)  # lewa kolumna – stała szerokość
        self.tab_opt.columnconfigure(1, weight=1)  # prawa kolumna – rośnie
        self.tab_opt.rowconfigure(0, weight=1)

        # === LEWA KOLUMNA: PARAMETRY ===
        opt_frame = ttk.LabelFrame(
            self.tab_opt,
            text="Parametry optymalizacji",
            padding=10
        )
        opt_frame.grid(row=0, column=0, sticky="ns", padx=10, pady=10)

        self.cso_vars = {
            "Karaluchy": tk.IntVar(value=30),
            "Iteracje": tk.IntVar(value=150),
            "Krok": tk.DoubleVar(value=0.5),
            "Max m": tk.IntVar(value=15),
            "C1 (koszt na jednostkę w systemie)": tk.DoubleVar(value=20.0),
            "C2 (koszt jednego kanału)": tk.DoubleVar(value=50.0),
        }

        for label, var in self.cso_vars.items():
            row = ttk.Frame(opt_frame)
            row.pack(fill="x", pady=4)
            ttk.Label(row, text=label, width=28, anchor="w").pack(side="left")
            ttk.Entry(row, textvariable=var, width=12).pack(side="left")

        ttk.Button(
            opt_frame,
            text="Optymalizuj (min. koszt całkowity)",
            command=self.run_opt
        ).pack(pady=15)

        # === PRAWA KOLUMNA ===
        right_frame = ttk.Frame(self.tab_opt)
        right_frame.grid(row=0, column=1, sticky="nsew", padx=10, pady=10)

        right_frame.columnconfigure(0, weight=1)
        right_frame.rowconfigure(0, weight=1)  # tabela
        right_frame.rowconfigure(1, weight=0)  # funkcja celu
        right_frame.rowconfigure(2, weight=0)  # efektywność / KPI

        # --- WIERSZ 1: TABELA ---
        self.table_opt = ttk.Treeview(
            right_frame,
            columns=("param", "val"),
            show="headings"
        )
        self.table_opt.heading("param", text="Parametr")
        self.table_opt.heading("val", text="Wartość")
        self.table_opt.column("param", width=450, anchor="w")
        self.table_opt.column("val", width=180, anchor="center")

        self.table_opt.grid(row=0, column=0, sticky="nsew")

        # --- WIERSZ 2: FUNKCJA CELU / KOSZT ---
        objective_frame = ttk.LabelFrame(right_frame, text="Funkcja celu", padding=8)
        objective_frame.grid(row=1, column=0, sticky="ew", pady=(10, 6))
        objective_frame.columnconfigure(0, weight=1)

        ttk.Label(
            objective_frame,
            text="J(m) = C1 · L(m) + C2 · m  (minimalizacja)",
            font=("Segoe UI", 11, "bold"),
            anchor="center"
        ).grid(row=0, column=0, sticky="ew")

        self.objective_label = tk.StringVar(
            value="J przed: -   |   J po: -   |   Oszczędność: -"
        )
        ttk.Label(
            objective_frame,
            textvariable=self.objective_label,
            font=("Segoe UI", 11),
            anchor="center"
        ).grid(row=1, column=0, sticky="ew", pady=(4, 0))

        self.cost_label = tk.StringVar(
            value="Optymalna liczba kanałów: -   |   Minimalny koszt całkowity: -"
        )
        ttk.Label(
            objective_frame,
            textvariable=self.cost_label,
            font=("Segoe UI", 11, "bold"),
            foreground="darkgreen",
            anchor="center"
        ).grid(row=2, column=0, sticky="ew", pady=(6, 0))

        # --- WIERSZ 3: EFEKTYWNOŚĆ / KPI ---
        kpi_frame = ttk.LabelFrame(right_frame, text="Efektywność (KPI)", padding=8)
        kpi_frame.grid(row=2, column=0, sticky="ew", pady=(0, 6))
        kpi_frame.columnconfigure(0, weight=1)

        self.efficiency_label = tk.StringVar(
            value="ρ: - → -   |   Wq: - → -   |   Lq: - → -"
        )
        ttk.Label(
            kpi_frame,
            textvariable=self.efficiency_label,
            font=("Segoe UI", 11),
            anchor="center"
        ).grid(row=0, column=0, sticky="ew")

    def run_opt(self):
        if not self.results_before:
            messagebox.showwarning("Uwaga", "Najpierw wykonaj analizę")
            return

        lam = self.lam_var.get()
        mu = self.mu_var.get()
        c1 = self.cso_vars["C1 (koszt na jednostkę w systemie)"].get()
        c2 = self.cso_vars["C2 (koszt jednego kanału)"].get()

        m_opt, total_cost = cockroach_optimization(
            lam, mu,
            self.cso_vars["Max m"].get(),
            c1, c2,
            self.cso_vars["Karaluchy"].get(),
            self.cso_vars["Iteracje"].get(),
            self.cso_vars["Krok"].get()
        )

        max_p = self.max_p_var.get()
        self.results_after = mm_m_full(lam, mu, m_opt, max_p)

        self.table_opt.delete(*self.table_opt.get_children())
        for k, v in self.results_after:
            if k == "" and v == "":
                self.table_opt.insert("", "end", values=("", ""), tags=("sep",))
            else:
                full_param = get_full_name(k)
                val_str = f"{v:.10f}" if isinstance(v, float) and not np.isinf(v) else "∞"
                self.table_opt.insert("", "end", values=(full_param, val_str))
        self.table_opt.tag_configure("sep", background="#ebebeb")
        # --- aktualizacja: funkcja celu + KPI (przed/po) ---
        before = self._results_to_dict(self.results_before)
        after = self._results_to_dict(self.results_after)

        m_before = int(round(float(before.get("m", 0))))
        m_after = int(round(float(after.get("m", 0))))

        L_before = float(before.get("L", 0))
        L_after = float(after.get("L", 0))

        J_before = c1 * L_before + c2 * m_before
        J_after = c1 * L_after + c2 * m_after
        saving = ((J_before - J_after) / J_before * 100) if J_before > 0 else 0

        rho_b = before.get("ρ", None)
        rho_a = after.get("ρ", None)
        Wq_b = before.get("Wq", None)
        Wq_a = after.get("Wq", None)
        Lq_b = before.get("Lq", None)
        Lq_a = after.get("Lq", None)

        self.objective_label.set(
            f"J przed: {self._fmt(J_before, 4)}   |   J po: {self._fmt(J_after, 4)}   |   Oszczędność: {self._fmt(saving, 2)}%"
        )

        self.efficiency_label.set(
            f"ρ: {self._fmt(rho_b, 4)} → {self._fmt(rho_a, 4)}   |   Wq: {self._fmt(Wq_b, 6)} → {self._fmt(Wq_a, 6)}   |   Lq: {self._fmt(Lq_b, 4)} → {self._fmt(Lq_a, 4)}"
        )

        self.cost_label.set(
            f"Optymalna liczba kanałów: {m_opt}   |   Minimalny koszt całkowity (C1·L + C2·m) = {total_cost:.4f}"
        )

    def compare_tab(self):
        # Układ siatki: tabela (góra) + podsumowanie (dół)
        self.tab_compare.columnconfigure(0, weight=1)
        self.tab_compare.rowconfigure(0, weight=1)  # tabela
        self.tab_compare.rowconfigure(1, weight=0)  # podsumowanie

        # --- TABELA PORÓWNAWCZA (zmniejszona wysokość) ---
        self.table_compare = ttk.Treeview(
            self.tab_compare,
            columns=("param", "before", "after"),
            show="headings",
            height=20
        )
        self.table_compare.heading("param", text="Parametr")
        self.table_compare.heading("before", text="Przed optymalizacją")
        self.table_compare.heading("after", text="Po optymalizacji")
        self.table_compare.column("param", width=450, anchor="w")
        self.table_compare.column("before", width=180, anchor="center")
        self.table_compare.column("after", width=180, anchor="center")
        self.table_compare.grid(row=0, column=0, sticky="nsew", padx=10, pady=(10, 6))

        # --- PODSUMOWANIE / KPI ---
        self.summary = tk.StringVar(value="")

        summary_frame = ttk.Frame(self.tab_compare)
        summary_frame.grid(row=1, column=0, sticky="ew", padx=10, pady=(0, 10))
        summary_frame.columnconfigure(0, weight=1)

        # wraplength: dopasowane do okna, żeby tekst nie uciekał poza widok
        self.summary_label = ttk.Label(
            summary_frame,
            textvariable=self.summary,
            font=("Segoe UI", 11, "bold"),
            anchor="w",
            justify="left",
            wraplength=1040
        )
        self.summary_label.grid(row=0, column=0, sticky="ew")

        # Aktualizuj automatycznie po wejściu na zakładkę
        self.tab_compare.bind("<Visibility>", lambda e: self.update_compare())

    def update_compare(self):
        if not self.results_before or not self.results_after:
            return

        self.table_compare.delete(*self.table_compare.get_children())

        for (k1, v1), (k2, v2) in zip(self.results_before, self.results_after):
            if k1 == "" and v1 == "":
                self.table_compare.insert("", "end", values=("", "", ""), tags=("sep",))
            else:
                full_param = get_full_name(k1)
                v1_str = f"{v1:.10f}" if isinstance(v1, float) and not np.isinf(v1) else "∞"
                v2_str = f"{v2:.10f}" if isinstance(v2, float) and not np.isinf(v2) else "∞"
                self.table_compare.insert("", "end", values=(full_param, v1_str, v2_str))

        self.table_compare.tag_configure("sep", background="#ebebeb")

        # Obliczenie kosztów przed i po
        c1 = self.cso_vars["C1 (koszt na jednostkę w systemie)"].get() if hasattr(self, 'cso_vars') else 0
        c2 = self.cso_vars["C2 (koszt jednego kanału)"].get() if hasattr(self, 'cso_vars') else 0

        L_before = [v for k, v in self.results_before if k == "L"][0]
        m_before = [v for k, v in self.results_before if k == "m"][0]
        cost_before = c1 * L_before + c2 * m_before

        L_after = [v for k, v in self.results_after if k == "L"][0]
        m_after = [v for k, v in self.results_after if k == "m"][0]
        cost_after = c1 * L_after + c2 * m_after

        saving = ((cost_before - cost_after) / cost_before * 100) if cost_before > 0 else 0


        # KPI efektywności (przed/po)
        before = self._results_to_dict(self.results_before)
        after = self._results_to_dict(self.results_after)

        rho_b = before.get("ρ", None)
        rho_a = after.get("ρ", None)
        Wq_b = before.get("Wq", None)
        Wq_a = after.get("Wq", None)
        Lq_b = before.get("Lq", None)
        Lq_a = after.get("Lq", None)

        self.summary.set(
            "Funkcja celu: J(m) = C1·L + C2·m\n"
            f"J przed: {cost_before:.4f} → J po: {cost_after:.4f}   |   Oszczędność: {saving:.2f}%   |   Kanały: {int(m_before)} → {int(m_after)}\n"
            f"Efektywność:  ρ: {self._fmt(rho_b, 4)} → {self._fmt(rho_a, 4)}   |   Wq: {self._fmt(Wq_b, 6)} → {self._fmt(Wq_a, 6)}   |   Lq: {self._fmt(Lq_b, 4)} → {self._fmt(Lq_a, 4)}"
        )

if __name__ == "__main__":
    app = QueueApp()
    app.mainloop()