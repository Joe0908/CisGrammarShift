#data resource: https://jaspar.elixir.no/temp/20251124234327_JASPAR2024_combined_matrices_1553058_jaspar.txt
import numpy as np
import matplotlib.pyplot as plt
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
import random
import operator 

# Setup and Seeding
# Set seed for reproducibility
np.random.seed(42)
random.seed(42)
torch.manual_seed(42)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(42)

# 1. DATASET CREATION (Probabilistic Biology)
def get_ctcf_probabilities():
    """
    Returns the Position Weight Matrix (PWM) for CTCF (MA0139.1).
    """
    pwm = np.array([
        [0.09, 0.32, 0.08, 0.51], [0.21, 0.18, 0.51, 0.10], [0.31, 0.06, 0.49, 0.14],
        [0.06, 0.88, 0.02, 0.04], [0.01, 0.99, 0.00, 0.00], [0.82, 0.01, 0.07, 0.10], 
        [0.04, 0.58, 0.37, 0.01], [0.12, 0.47, 0.05, 0.36], [0.93, 0.01, 0.04, 0.02], 
        [0.01, 0.00, 0.99, 0.00], [0.37, 0.00, 0.62, 0.01], [0.06, 0.01, 0.55, 0.38], 
        [0.01, 0.01, 0.98, 0.00], [0.06, 0.00, 0.86, 0.08], [0.11, 0.81, 0.01, 0.07],
        [0.41, 0.01, 0.56, 0.02], [0.09, 0.54, 0.34, 0.03], [0.13, 0.36, 0.08, 0.43],
        [0.44, 0.20, 0.29, 0.07]
    ])
    return pwm

def generate_probabilistic_data(num_sequences=4000, seq_length=100):
    sequences = []
    labels = []
    bases = ['A', 'C', 'G', 'T']
    pwm = get_ctcf_probabilities()
    motif_len = len(pwm)
    
    print("Generating data using Probabilistic CTCF Matrix...")

    for i in range(num_sequences):
        seq_list = random.choices(bases, k=seq_length)
        label = 1 if i < num_sequences // 2 else 0
        
        if label == 1:
            motif_instance = [np.random.choice(bases, p=probs) for probs in pwm]
            insert_loc = random.randint(0, seq_length - motif_len)
            for j, char in enumerate(motif_instance):
                seq_list[insert_loc + j] = char

        sequences.append("".join(seq_list))
        labels.append(label)

    return np.array(sequences), np.array(labels)

# 2. ONE-HOT ENCODING
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

# 3. MODELS IN PYTORCH 
class CNN(nn.Module):
    def __init__(self):
        super().__init__()
        self.conv = nn.Conv1d(4, 32, kernel_size=12) 
        self.pool = nn.AdaptiveMaxPool1d(1)
        self.fc1 = nn.Linear(32, 16)
        self.fc2 = nn.Linear(16, 1)

    def forward(self, x, return_logits=False): 
        x = x.permute(0, 2, 1) 
        x = torch.relu(self.conv(x))
        x = self.pool(x).squeeze(-1)
        x = torch.relu(self.fc1(x))
        logits = self.fc2(x)
        if return_logits:
            return logits
        return torch.sigmoid(logits)

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
        h = torch.cat((h[-2], h[-1]), dim=1) 
        x = torch.relu(self.fc1(h))
        return torch.sigmoid(self.fc2(x))

class TransformerModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.embed = nn.Linear(4, 32)
        self.pos_embedding = nn.Parameter(torch.randn(1, SEQ_LEN, 32))
        
        encoder_layer = nn.TransformerEncoderLayer(d_model=32, nhead=2, batch_first=True)
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers=1)
        self.fc1 = nn.Linear(32, 16)
        self.fc2 = nn.Linear(16, 1)

    def forward(self, x):
        x = self.embed(x)
        x = x + self.pos_embedding 
        x = self.encoder(x)
        x = x.mean(dim=1) 
        x = torch.relu(self.fc1(x))
        return torch.sigmoid(self.fc2(x))

# 4. DATA PREPARATION (MANUAL SPLIT) 
raw_seqs, labels = generate_probabilistic_data(seq_length=SEQ_LEN)
X = one_hot_encode(raw_seqs)
y = labels.astype(np.float32)

# MANUAL SPLITTING LOGIC (No sklearn needed)
indices = np.arange(len(X))
np.random.shuffle(indices) 

train_split = int(0.7 * len(X))
val_split = int(0.85 * len(X))

train_idx = indices[:train_split]
val_idx = indices[train_split:val_split]
test_idx = indices[val_split:]

X_train, y_train = X[train_idx], y[train_idx]
X_val, y_val = X[val_idx], y[val_idx]
X_test, y_test = X[test_idx], y[test_idx]

print(f"Data shapes: Train {X_train.shape}, Val {X_val.shape}, Test {X_test.shape}")

# Loaders
train_loader = DataLoader(TensorDataset(torch.tensor(X_train), torch.tensor(y_train)), batch_size=32, shuffle=True)
val_loader = DataLoader(TensorDataset(torch.tensor(X_val), torch.tensor(y_val)), batch_size=32)
test_loader = DataLoader(TensorDataset(torch.tensor(X_test), torch.tensor(y_test)), batch_size=32)

# 5. TRAINING FUNCTION
def train_model(model):
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    loss_fn = nn.BCELoss()
    history = {"val_acc": [], "val_loss": []}

    print(f"Training {model.__class__.__name__}...")
    for epoch in range(15):
        model.train()
        for xb, yb in train_loader:
            if isinstance(model, CNN):
                preds = model(xb).squeeze()
            else:
                preds = model(xb).squeeze()
            
            loss = loss_fn(preds, yb)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        # Validation
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

        val_acc = correct / total
        history["val_loss"].append(val_loss / len(val_loader))
        history["val_acc"].append(val_acc)
        
        if (epoch + 1) % 5 == 0:
            print(f"  Epoch {epoch+1}: Val Acc = {val_acc:.4f}")

    return history

# Saliency Map Visualisation Function 
def visualize_saliency_map(model, one_hot_sequence, true_label, title="Saliency Map"):
    # NOTE: This function only works correctly for the CNN model structure
    model.eval()
    
    input_tensor = torch.tensor(one_hot_sequence, dtype=torch.float32).unsqueeze(0)
    input_tensor.requires_grad_(True)
    
    output = model(input_tensor, return_logits=True) 
    
    output.backward(torch.ones_like(output))
    
    gradients = input_tensor.grad.data.squeeze(0).cpu().numpy()
    
    # Saliency = Sum of absolute gradients across nucleotide channels
    saliency_scores = np.sum(np.abs(gradients), axis=1)
    
    # Normalise scores for visualisation
    saliency_scores = (saliency_scores - np.min(saliency_scores)) / (np.max(saliency_scores) - np.min(saliency_scores) + 1e-8)
    
    # Convert one-hot sequence back to string
    bases = ['A', 'C', 'G', 'T']
    original_sequence = "".join([bases[np.argmax(pos_vector)] for pos_vector in one_hot_sequence])

    # Plotting
    plt.figure(figsize=(15, 3))
    plt.imshow(saliency_scores[np.newaxis, :], cmap='hot', aspect='auto') 
    cbar = plt.colorbar(label='Saliency Score (Importance)')
    cbar.set_ticks([0, 0.5, 1])
    plt.yticks([]) 
    plt.xticks(np.arange(0, SEQ_LEN, 5)) 
    
    # Annotate sequence
    for i, char in enumerate(original_sequence):
        plt.text(i, -0.4, char, ha='center', va='center', color='black', fontsize=8) 

    plt.title(f"{title} (True Label: {int(true_label)})")
    plt.xlabel("Sequence Position")
    plt.tight_layout()
    plt.show()

#  6. RUN MODELS 
models = {
    "CNN": CNN(),
    "SimpleRNN": SimpleRNN(),
    "LSTM": LSTMModel(),
    "Transformer": TransformerModel()
}

histories = {}
results = {}

for name, model in models.items():
    hist = train_model(model)
    histories[name] = hist
    
    # Test Evaluation
    model.eval()
    correct = 0
    total = 0
    with torch.no_grad():
        for xb, yb in test_loader:
            # Note: For Saliency Map compatibility, ensure the output here is sigmoid activated
            preds = (model(xb).squeeze() > 0.5).float() 
            correct += (preds == yb).sum().item()
            total += len(yb)
    results[name] = correct / total

# 7. RANKING AND PLOTS 

# Generate Line Charts for Validation Metrics
try:
    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    for name, hist in histories.items():
        plt.plot(hist["val_acc"], label=name)
    plt.title("Validation Accuracy")
    plt.xlabel("Epochs")
    plt.legend()

    plt.subplot(1, 2, 2)
    for name, hist in histories.items():
        plt.plot(hist["val_loss"], label=name)
    plt.title("Validation Loss")
    plt.xlabel("Epochs")
    plt.legend()
    plt.tight_layout()
    plt.show()
except ImportError:
    print("Matplotlib not found, skipping plots.")

# Rank Models based on Test Accuracy
ranked_results = sorted(results.items(), key=operator.itemgetter(1), reverse=True)

print("\n--- FINAL MODEL RANKING (Based on Test Accuracy) ---")
print(f"{'Rank':<5} | {'Model':<15} | {'Accuracy':<10}")
print("-" * 35)
for i, (name, acc) in enumerate(ranked_results):
    print(f"{i+1:<5} | {name:<15} | {acc:.4f}")

# 8. Saliency Map Visualisation (Dynamically selects Rank 1 model) 

best_model_name = ranked_results[0][0]
best_model = models[best_model_name]

if best_model_name == "CNN":
    print(f"\n--- Generating Saliency Maps for Best Model ({best_model_name}) ---")

    # Find a correctly classified positive and negative example from the test set
    # Using the best model for prediction
    best_model_preds = (best_model(torch.tensor(X_test, dtype=torch.float32)).squeeze() > 0.5).float().cpu().numpy()
    correct_pos_idx = np.where((y_test == 1) & (best_model_preds == 1))[0][0] 
    correct_neg_idx = np.where((y_test == 0) & (best_model_preds == 0))[0][0] 

    pos_seq_onehot = X_test[correct_pos_idx]
    neg_seq_onehot = X_test[correct_neg_idx]
    pos_true_label = y_test[correct_pos_idx]
    neg_true_label = y_test[correct_neg_idx]

    # Visualise for a correctly classified positive example (Motif Present)
    visualize_saliency_map(best_model, pos_seq_onehot, pos_true_label, 
                           title=f"{best_model_name} Saliency Map - Positive Example (Motif Present) ")

    # Visualise for a correctly classified negative example (Random Noise)
    visualize_saliency_map(best_model, neg_seq_onehot, neg_true_label, 
                           title=f"{best_model_name} Saliency Map - Negative Example (No Motif) ")
else:
    print(f"\nSaliency Map skipped: Best model is the {best_model_name}, but visualisation is only supported for the CNN structure.")
