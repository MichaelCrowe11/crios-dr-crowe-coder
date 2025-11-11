# Vercel Serverless Functions (Python)

This directory contains Python serverless functions for Vercel deployment.

## ⚠️ Important Limitation

**RDKit is too large for Vercel serverless functions** (~300MB, Vercel limit is 250MB).

## Recommended Solutions

### Option 1: Hybrid Architecture (Recommended)
Deploy the frontend on Vercel and keep the FastAPI backend with RDKit on Railway:

```
Frontend (Vercel) → Backend API (Railway)
                    └─ RDKit, FastAPI, Celery
```

**Pros:**
- Full RDKit functionality
- Fast frontend deployment
- Best of both platforms

**Setup:**
1. Deploy frontend to Vercel
2. Deploy backend to Railway (see `RAILWAY_DEPLOYMENT.md`)
3. Set `NEXT_PUBLIC_API_URL` to Railway backend URL

### Option 2: External Chemistry APIs
Use public chemistry APIs instead of RDKit:

**PubChem API:**
```python
# Validate SMILES via PubChem
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/JSON
```

**ChemSpider:**
```python
# Requires API key from https://developer.rsc.org/
POST https://api.rsc.org/compounds/v1/filter/smiles
```

**CDK (Chemistry Development Kit):**
- Java-based chemistry library
- Can be wrapped in Node.js serverless function
- Smaller than RDKit

### Option 3: Lightweight Validation
Use the lightweight JavaScript/Python validators included in this project:

- `/pages/api/validate.ts` - TypeScript validation (basic)
- `/api/chemistry/validate.py` - Python validation (basic)

**Note:** These are NOT suitable for production chemistry applications.

## Current Implementation

The Python functions in this directory use lightweight validation as a **proof of concept**.

For production use, implement Option 1 (Hybrid Architecture).

## File Structure

```
api/
├── README.md                 # This file
├── requirements.txt          # Python dependencies (minimal)
└── chemistry/
    └── validate.py          # Lightweight validation example
```

## Testing Locally

```bash
# Install Vercel CLI
npm i -g vercel

# Run local development server
vercel dev

# Test Python function
curl -X POST http://localhost:3000/api/chemistry/validate \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CCO"}'
```

## Deployment

These functions are automatically deployed when you deploy to Vercel.

```bash
vercel --prod
```

## Recommendations

For the CriOS Discovery Engine, we recommend:

1. **Frontend on Vercel** - Fast, global CDN, excellent Next.js support
2. **Backend on Railway** - Full RDKit, no size limits, persistent services
3. **Database on Vercel Postgres or Neon** - Serverless PostgreSQL
4. **Cache on Upstash Redis** - Serverless Redis

This gives you the best performance and functionality.

---

© 2025 CriOS Discovery Engine
