# CriOS Nova Dashboard UI

A modern React/Next.js frontend for the CriOS Discovery Engine with live API integration.

## Features

- **Real-time API Integration**: Connects to CriOS backend for live molecular validation, scoring, and ethics checks
- **Interactive Compound Discovery**: Advanced filtering, Pareto optimization, and candidate prioritization
- **Crowe Discovery Framework**: Live scoring with multi-objective optimization
- **Ethics Compliance**: Real-time compliance checking against ethical policies
- **3D Visualization**: Ready for molecular structure and protein visualization
- **Performance Optimized**: Virtual scrolling, web workers, and debounced operations

## Quick Start

```bash
# Install dependencies
npm install

# Start development server
npm run dev

# Open http://localhost:3000
```

## API Integration

The dashboard connects to the CriOS backend API at `http://localhost:8000` by default.

### Available Endpoints:
- `POST /validate` - Molecular validation and standardization
- `POST /crowe-score` - Crowe Discovery Framework scoring
- `POST /ethics-check` - Ethics compliance checking  
- `POST /similarity` - Molecular similarity search
- `GET /config` - System configuration
- `GET /health` - API health status

## Architecture

```
ui/
├── app/                 # Next.js 13+ App Router
├── components/ui/       # shadcn/ui components
├── nova-dashboard.tsx   # Main dashboard component
└── lib/                # Utility functions
```

## Technology Stack

- **Framework**: Next.js 14 with React 18
- **Styling**: Tailwind CSS with shadcn/ui components
- **Charts**: Recharts for data visualization
- **Animations**: Framer Motion
- **Icons**: Lucide React
- **TypeScript**: Full type safety

## Key Components

### NovaDashboard
Main dashboard component with live API integration, real-time updates, and comprehensive molecular analysis.

### API Service
Centralized API client with error handling and type safety:
```typescript
const apiService = {
  validateMolecule: async (smiles: string) => { /* ... */ },
  croweScore: async (smiles: string) => { /* ... */ },
  ethicsCheck: async (smiles: string) => { /* ... */ }
};
```

### Performance Features
- Virtual scrolling for large datasets
- Web Workers for Pareto computations
- Debounced filters and search
- Optimistic UI updates

## Development

```bash
# Development server
npm run dev

# Type checking
npm run type-check

# Build for production
npm run build

# Start production server
npm start
```

## Environment Variables

```bash
CRIOS_API_URL=http://localhost:8000  # CriOS backend URL
```

## Integration with CriOS Backend

1. **Start CriOS API**: `python -m uvicorn crios.web.main:app --reload --port 8000`
2. **Start Frontend**: `npm run dev`
3. **Access Dashboard**: `http://localhost:3000`

The UI will automatically connect to the CriOS API and provide real-time molecular analysis capabilities.

## Features in Development

- [ ] 3D molecular structure visualization (3Dmol.js/NGL)
- [ ] Protein-ligand interaction viewer
- [ ] Advanced synthesis route visualization
- [ ] Real-time collaborative features
- [ ] Export to laboratory management systems

---

**CriOS Discovery Engine** • Science before status. Discovery before profit. • © 2025 Crowe BioSystems