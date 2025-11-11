# CriOS Discovery Engine - Vercel Deployment Guide

Complete guide for deploying the CriOS Discovery Engine to Vercel.

## 🚀 Deployment Options

### Option 1: Hybrid Architecture (Recommended for Production)

Deploy frontend on Vercel, backend on Railway for full RDKit support.

```
┌─────────────────────────────────────────┐
│         Vercel (Frontend)                │
│  - Next.js UI (Global CDN)              │
│  - Static pages                          │
│  - Client-side React                     │
└─────────────┬───────────────────────────┘
              │ API calls via HTTPS
              ▼
┌─────────────────────────────────────────┐
│        Railway (Backend)                 │
│  - FastAPI + RDKit                      │
│  - PostgreSQL Database                   │
│  - Redis Cache                           │
│  - Celery Workers                        │
└─────────────────────────────────────────┘
```

**Pros:**
- ✅ Full RDKit functionality
- ✅ Fast global CDN for frontend
- ✅ No serverless size limits
- ✅ Background job support (Celery)
- ✅ Best performance

**Cost:** $5-10/month (Vercel) + $25-45/month (Railway) = **$30-55/month**

### Option 2: Pure Vercel (Limited Functionality)

Deploy everything on Vercel with lightweight chemistry validation.

```
┌─────────────────────────────────────────┐
│            Vercel Only                   │
│  - Next.js Frontend                     │
│  - Next.js API Routes (lightweight)     │
│  - External DB (Vercel Postgres)        │
│  - External Cache (Upstash Redis)       │
└─────────────────────────────────────────┘
```

**Limitations:**
- ❌ No RDKit (too large for serverless)
- ❌ No Celery workers (long-running processes)
- ⚠️ Basic SMILES validation only
- ⚠️ External API calls needed for chemistry

**Cost:** $20-30/month (hobby tier)

---

## 📋 Option 1: Hybrid Deployment (Recommended)

### Step 1: Deploy Backend to Railway

Follow the [Railway Deployment Guide](RAILWAY_DEPLOYMENT.md) to deploy:
- FastAPI backend with RDKit
- PostgreSQL database
- Redis cache
- Celery workers

```bash
# Quick Railway deployment
./deploy/railway-deploy.sh
```

### Step 2: Deploy Frontend to Vercel

#### A. Using Vercel CLI

```bash
# Install Vercel CLI
npm install -g vercel

# Navigate to frontend directory
cd platform/frontend

# Login to Vercel
vercel login

# Deploy
vercel --prod
```

#### B. Using Vercel Dashboard

1. Go to [Vercel Dashboard](https://vercel.com/new)
2. Click **"Import Project"**
3. Select your GitHub repository
4. Configure project:
   - **Framework Preset**: Next.js
   - **Root Directory**: `platform/frontend`
   - **Build Command**: `npm run build`
   - **Output Directory**: `.next`

5. Set environment variables:

```bash
# Frontend environment variables
NEXT_PUBLIC_API_URL=https://your-backend.railway.app
NEXT_PUBLIC_WS_URL=wss://your-backend.railway.app
NODE_ENV=production
```

6. Click **"Deploy"**

### Step 3: Configure Backend CORS

Update your Railway backend to allow Vercel frontend:

```bash
# In Railway dashboard, add environment variable:
CORS_ORIGINS=https://your-project.vercel.app
```

Redeploy the backend service.

### Step 4: Verify Deployment

```bash
# Check frontend
curl https://your-project.vercel.app

# Check API connection
curl https://your-project.vercel.app/api/health

# Should return backend health status via proxy
```

---

## 📋 Option 2: Pure Vercel Deployment

### Prerequisites

- Vercel account
- GitHub repository
- Upstash Redis account (free tier available)
- Vercel Postgres (optional, or use Neon/Supabase)

### Step 1: Set Up External Services

#### A. Upstash Redis

1. Go to [Upstash Console](https://console.upstash.com/)
2. Create new Redis database
3. Copy connection details:
   - `REDIS_URL`
   - `UPSTASH_REDIS_REST_URL`
   - `UPSTASH_REDIS_REST_TOKEN`

#### B. Database (Choose One)

**Option A: Vercel Postgres**
```bash
# In Vercel dashboard
# Storage → Create Database → Postgres
# Copy environment variables automatically
```

**Option B: Neon (Serverless Postgres)**
1. Go to [Neon Console](https://console.neon.tech/)
2. Create new project
3. Copy `DATABASE_URL`

**Option C: Supabase**
1. Go to [Supabase](https://supabase.com/)
2. Create new project
3. Copy connection string from Settings → Database

### Step 2: Deploy to Vercel

#### Using Vercel CLI

```bash
# From repository root
vercel

# Follow prompts:
# - Link to existing project or create new
# - Set root directory: platform/frontend
# - Override settings: Yes
# - Build Command: npm run build
# - Output Directory: .next
# - Development Command: npm run dev

# Add environment variables
vercel env add NEXT_PUBLIC_API_URL
# Enter: https://your-project.vercel.app

vercel env add REDIS_URL
# Enter your Upstash Redis URL

vercel env add POSTGRES_URL
# Enter your database URL

vercel env add ANTHROPIC_API_KEY
# Enter your Anthropic API key

# Deploy to production
vercel --prod
```

#### Using GitHub Integration

1. Go to [Vercel Dashboard](https://vercel.com/new)
2. Import your GitHub repository
3. Configure settings:

```
Framework Preset: Next.js
Root Directory: platform/frontend
Build Command: npm run build
Output Directory: .next
Install Command: npm install
```

4. Add environment variables (see Environment Variables section below)
5. Click **"Deploy"**

### Step 3: Configure Environment Variables

In Vercel Dashboard → Settings → Environment Variables:

```bash
# Public variables (embedded in build)
NEXT_PUBLIC_API_URL=https://your-project.vercel.app
NEXT_PUBLIC_WS_URL=wss://your-project.vercel.app

# Server-side variables
REDIS_URL=redis://default:***@***-us1-alive-mallard-12345.upstash.io:6379
UPSTASH_REDIS_REST_URL=https://***-us1-alive-mallard-12345.upstash.io
UPSTASH_REDIS_REST_TOKEN=***

# Database (if using Vercel Postgres, these are auto-set)
POSTGRES_URL=postgres://***
POSTGRES_PRISMA_URL=postgres://***?pgbouncer=true
POSTGRES_URL_NON_POOLING=postgres://***

# API Keys
ANTHROPIC_API_KEY=sk-ant-***

# Security
SECRET_KEY=<generate-with-openssl-rand-base64-32>

# Optional: Background jobs with QStash
QSTASH_URL=https://qstash.upstash.io
QSTASH_TOKEN=***
```

### Step 4: Update Configuration

The deployment uses `next.config.vercel.js`. Copy it to production:

```bash
cd platform/frontend
cp next.config.vercel.js next.config.js
git add next.config.js
git commit -m "Update Next.js config for Vercel"
git push
```

Vercel will auto-deploy on push.

### Step 5: Add Integrations

#### Vercel Postgres Integration

```bash
# From Vercel dashboard
# Storage → Add → Postgres
# This automatically sets environment variables
```

#### Upstash Redis Integration

```bash
# From Vercel dashboard
# Integrations → Browse Marketplace → Upstash Redis
# Connect your Upstash account
# Select Redis database
# Auto-configures environment variables
```

---

## 🔧 Project Structure for Vercel

```
crios-dr-crowe-coder/
├── platform/frontend/          # Next.js frontend
│   ├── pages/
│   │   ├── index.tsx          # Landing page
│   │   ├── ide.tsx            # IDE interface
│   │   ├── pricing.tsx        # Pricing page
│   │   └── api/               # Next.js API routes
│   │       ├── health.ts      # Health check
│   │       ├── docs.ts        # API documentation
│   │       └── validate.ts    # Basic validation
│   ├── components/            # React components
│   ├── package.json
│   ├── next.config.js         # Next.js configuration
│   └── tsconfig.json
├── api/                       # Python serverless functions
│   ├── chemistry/
│   │   └── validate.py        # Lightweight validation
│   ├── requirements.txt       # Python deps (minimal)
│   └── README.md
├── vercel.json                # Vercel project config
└── .env.vercel.example        # Environment template
```

---

## 📊 API Routes Architecture

### Next.js API Routes (TypeScript)

Located in `platform/frontend/pages/api/`:

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/health` | GET | Health check |
| `/api/docs` | GET | API documentation |
| `/api/validate` | POST | Basic SMILES validation |

### Python Serverless Functions

Located in `/api/`:

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/chemistry/validate` | POST | Python-based validation (lightweight) |

### For Full Chemistry Features

Use hybrid architecture with Railway backend for:
- RDKit molecule processing
- Similarity search
- Crowe scoring
- Complex computational chemistry

---

## 🔒 Security Configuration

### Environment Variables

Never commit secrets to git. Use Vercel's encrypted environment storage.

### CORS Configuration

For hybrid deployment, configure CORS in backend:

```python
# platform/backend/app.py
from fastapi.middleware.cors import CORSMiddleware

app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "https://your-project.vercel.app",
        "https://www.your-domain.com"
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
```

### API Keys

Store in Vercel environment variables:
- `ANTHROPIC_API_KEY`
- Database credentials
- Redis credentials

---

## 📈 Performance Optimization

### Edge Functions

For ultra-low latency, convert API routes to Edge Functions:

```typescript
// pages/api/edge-example.ts
export const config = {
  runtime: 'edge',
};

export default async function handler(req: Request) {
  return new Response('Hello from the Edge!');
}
```

### Caching Strategy

```typescript
// next.config.js
module.exports = {
  async headers() {
    return [
      {
        source: '/api/:path*',
        headers: [
          {
            key: 'Cache-Control',
            value: 's-maxage=60, stale-while-revalidate',
          },
        ],
      },
    ];
  },
};
```

### Image Optimization

```typescript
// next.config.js
module.exports = {
  images: {
    domains: ['your-backend.railway.app'],
    deviceSizes: [640, 750, 828, 1080, 1200, 1920, 2048, 3840],
    imageSizes: [16, 32, 48, 64, 96, 128, 256, 384],
    formats: ['image/webp'],
  },
};
```

---

## 🛠️ Troubleshooting

### Build Failures

**Error: Module not found**
```bash
# Check package.json in platform/frontend
cd platform/frontend
npm install
```

**Error: Environment variable not set**
```bash
# Verify in Vercel Dashboard → Settings → Environment Variables
# Or use Vercel CLI:
vercel env pull .env.local
```

### Runtime Errors

**API calls failing**
```bash
# Check NEXT_PUBLIC_API_URL is set correctly
# Verify CORS is configured on backend
# Check backend health: curl https://backend-url/health
```

**Database connection issues**
```bash
# Verify DATABASE_URL environment variable
# Check connection string format
# Test connection from Vercel Functions:
vercel logs
```

### Performance Issues

**Slow cold starts**
```bash
# Use Edge Functions for critical paths
# Enable caching headers
# Optimize bundle size: npm run analyze
```

**High memory usage**
```bash
# Increase function memory in vercel.json
# Optimize React rendering
# Use dynamic imports
```

---

## 📊 Monitoring & Analytics

### Vercel Analytics

Enable in Vercel Dashboard:
```bash
# Settings → Analytics → Enable
```

Add to `_app.tsx`:
```typescript
import { Analytics } from '@vercel/analytics/react';

export default function App({ Component, pageProps }) {
  return (
    <>
      <Component {...pageProps} />
      <Analytics />
    </>
  );
}
```

### Logging

View logs in real-time:
```bash
vercel logs --follow

# Or specific deployment:
vercel logs [deployment-url]
```

### Error Tracking

Integrate Sentry:
```bash
npm install --save @sentry/nextjs

# Run setup wizard:
npx @sentry/wizard -i nextjs
```

---

## 💰 Cost Estimation

### Hobby Plan (Free Tier Included)

| Service | Cost |
|---------|------|
| Vercel Hobby | $0 (100GB bandwidth, 100 deployments/day) |
| Upstash Redis | $0 (10K commands/day) |
| Neon Postgres | $0 (0.5GB storage, 100 hours compute) |
| **Total** | **$0-10/month** |

### Pro Plan (Production)

| Service | Cost |
|---------|------|
| Vercel Pro | $20/month |
| Upstash Redis Pro | $10/month |
| Neon Scale | $19/month |
| Vercel Postgres | Included |
| **Total** | **$50-70/month** |

### Hybrid (Vercel + Railway)

| Service | Cost |
|---------|------|
| Vercel (Frontend) | $0-20/month |
| Railway (Backend) | $25-45/month |
| **Total** | **$25-65/month** |

---

## 🚦 Going to Production

### Pre-Launch Checklist

- [ ] Environment variables configured
- [ ] Custom domain added and verified
- [ ] SSL certificate active (automatic)
- [ ] Analytics enabled
- [ ] Error tracking configured
- [ ] Performance monitoring set up
- [ ] Database backups enabled
- [ ] CORS properly configured
- [ ] API rate limiting implemented
- [ ] Security headers configured

### Custom Domain

```bash
# Add domain in Vercel Dashboard:
# Settings → Domains → Add Domain

# Add DNS records:
# Type: CNAME
# Name: www or @
# Value: cname.vercel-dns.com
```

### SSL/HTTPS

Vercel automatically provisions SSL certificates via Let's Encrypt.

---

## 📚 Additional Resources

- **Vercel Documentation**: https://vercel.com/docs
- **Next.js Documentation**: https://nextjs.org/docs
- **Vercel Templates**: https://vercel.com/templates
- **Railway Deployment**: See [RAILWAY_DEPLOYMENT.md](RAILWAY_DEPLOYMENT.md)
- **CriOS Documentation**: See [README.md](README.md)

---

## 🎉 Quick Start Commands

```bash
# Install Vercel CLI
npm i -g vercel

# Login
vercel login

# Deploy to production
cd platform/frontend
vercel --prod

# View logs
vercel logs --follow

# Open dashboard
vercel open

# Pull environment variables
vercel env pull

# Add domain
vercel domains add your-domain.com
```

---

## 💡 Recommendations

For the **CriOS Discovery Engine**, we recommend:

### Development
- Local development with `npm run dev`
- Use Vercel CLI for testing: `vercel dev`

### Staging
- Deploy to Vercel preview deployments (automatic on PR)
- Test with production-like environment

### Production
- **Option 1 (Recommended)**: Hybrid Vercel (frontend) + Railway (backend)
- **Option 2 (Lightweight)**: Pure Vercel with external APIs for chemistry

This ensures:
- ✅ Fast global CDN for frontend
- ✅ Full RDKit functionality
- ✅ Scalable architecture
- ✅ Cost-effective deployment
- ✅ Easy maintenance

---

**CriOS Discovery Engine** • Deployed with Vercel • Science before status. Discovery before profit.

© 2025 Crowe BioSystems
