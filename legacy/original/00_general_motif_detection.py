import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from sklearn.model_selection import train_test_split
import random

# --- 1. DATASET CREATION ---
def generate_synthetic_data(num_sequences=5000, seq_length=100, motif="ATGCCA"):
    bases = ['A', 'C', 'G', 'T']
    sequences, labels = [], []

    for i in range(num_sequences):
        seq = "".join(random.choices(bases, k=seq_length))
        label = 1 if i < num_sequences // 2 else 0
        
        if label == 1:
            insert_loc = random.randint(0, seq_length - len(motif))
            seq = seq[:insert_loc] + motif + seq[insert_loc+len(motif):]

        sequences.append(seq)
        labels.append(label)

    combined = list(zip(sequences, labels))
    random.shuffle(combined)
    sequences, labels = zip(*combined)
    return np.array(sequences), np.array(labels)

# --- 2. ONE-HOT ENCODING ---
def one_hot_encode(sequences):
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    encoded_list = []

    for seq in sequences:
        mat = np.zeros((len(seq), 4), dtype=np.float32)
        for i, nt in enumerate(seq):
            mat[i, mapping[nt]] = 1
        encoded_list.append(mat)

    return np.array(encoded_list)

SEQ_LEN = 100
NUM_BASES = 4

# --- 3. MODELS IN PYTORCH ---
class CNN(nn.Module):
    def __init__(self):
        super().__init__()
        self.conv = nn.Conv1d(4, 32, kernel_size=12)
        self.pool = nn.AdaptiveMaxPool1d(1)
        self.fc1 = nn.Linear(32, 16)
        self.fc2 = nn.Linear(16, 1)

    def forward(self, x):
        x = x.permute(0, 2, 1)
        x = torch.relu(self.conv(x))
        x = self.pool(x).squeeze(-1)
        x = torch.relu(self.fc1(x))
        return torch.sigmoid(self.fc2(x))

class SimpleRNN(nn.Module):
    def __init__(self):
        super().__init__()
        self.rnn = nn.RNN(4, 32, batch_first=True)
        self.fc1 = nn.Linear(32, 16)
        self.fc2 = nn.Linear(16, 1)

    def forward(self, x):
        _, h = self.rnn(x)
        x = torch.relu(self.fc1(h.squeeze(0)))
        return torch.sigmoid(self.fc2(x))

class LSTMModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.lstm = nn.LSTM(4, 16, batch_first=True, bidirectional=True)
        self.fc1 = nn.Linear(32, 16)
        self.fc2 = nn.Linear(16, 1)

    def forward(self, x):
        _, (h, _) = self.lstm(x)
        h = torch.cat((h[0], h[1]), dim=1)
        x = torch.relu(self.fc1(h))
        return torch.sigmoid(self.fc2(x))

class TransformerModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.embed = nn.Linear(4, 32)
        encoder_layer = nn.TransformerEncoderLayer(d_model=32, nhead=2)
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers=1)
        self.fc1 = nn.Linear(32, 16)
        self.fc2 = nn.Linear(16, 1)

    def forward(self, x):
        x = self.embed(x)
        x = self.encoder(x)
        x = x.mean(dim=1)
        x = torch.relu(self.fc1(x))
        return torch.sigmoid(self.fc2(x))

# --- 4. DATA PREPARATION ---
raw_seqs, labels = generate_synthetic_data(num_sequences=4000, seq_length=SEQ_LEN)
X = one_hot_encode(raw_seqs).astype(np.float32)
y = labels.astype(np.float32)

X_train, X_tmp, y_train, y_tmp = train_test_split(X, y, test_size=0.3, random_state=42)
X_val, X_test, y_val, y_test = train_test_split(X_tmp, y_tmp, test_size=0.5, random_state=42)

train_loader = DataLoader(TensorDataset(torch.tensor(X_train), torch.tensor(y_train)), batch_size=32, shuffle=True)
val_loader = DataLoader(TensorDataset(torch.tensor(X_val), torch.tensor(y_val)), batch_size=32)
test_loader = DataLoader(TensorDataset(torch.tensor(X_test), torch.tensor(y_test)), batch_size=32)

# --- 5. TRAINING FUNCTION ---
def train_model(model):
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    loss_fn = nn.BCELoss()

    history = {"val_acc": [], "val_loss": []}

    for epoch in range(15):
        model.train()
        for xb, yb in train_loader:
            preds = model(xb).squeeze()
            loss = loss_fn(preds, yb)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        # validation
        model.eval()
        correct = 0
        total = 0
        val_loss = 0

        with torch.no_grad():
            for xb, yb in val_loader:
                preds = model(xb).squeeze()
                val_loss += loss_fn(preds, yb).item()

                preds_bin = (preds > 0.5).float()
                correct += (preds_bin == yb).sum().item()
                total += len(yb)

        history["val_loss"].append(val_loss / len(val_loader))
        history["val_acc"].append(correct / total)

        print(f"Epoch {epoch+1}: Val Acc = {correct/total:.4f}")

    return history

# --- 6. RUN MODELS ---
models = {
    "CNN": CNN(),
    "SimpleRNN": SimpleRNN(),
    "LSTM": LSTMModel(),
    "Transformer": TransformerModel()
}

histories = {}
results = {}

for name, model in models.items():
    print(f"\n--- Training {name} ---")
    hist = train_model(model)
    histories[name] = hist

    # test accuracy
    model.eval()
    correct = 0
    total = 0
    with torch.no_grad():
        for xb, yb in test_loader:
            preds = (model(xb).squeeze() > 0.5).float()
            correct += (preds == yb).sum().item()
            total += len(yb)
    results[name] = correct / total

print("\n--- Final Test Accuracy ---")
print(pd.DataFrame(results, index=["Accuracy"]).T)

# --- 7. PLOTS ---
plt.figure(figsize=(12, 5))
for name, hist in histories.items():
    plt.plot(hist["val_acc"], label=name)
plt.title("Validation Accuracy")
plt.legend()

plt.figure(figsize=(12, 5))
for name, hist in histories.items():
    plt.plot(hist["val_loss"], label=name)
plt.title("Validation Loss")
plt.legend()
plt.show()
