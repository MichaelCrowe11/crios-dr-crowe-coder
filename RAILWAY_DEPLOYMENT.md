# CriOS Discovery Engine - Railway Deployment Guide

Complete guide for deploying the CriOS Discovery Engine to Railway.

## 🚀 Quick Deploy (Recommended)

### Prerequisites
- Railway account (sign up at https://railway.app)
- GitHub repository connected to Railway
- Anthropic API key

### One-Click Deployment

[![Deploy on Railway](https://railway.app/button.svg)](https://railway.app/template)

Or deploy manually using the steps below.

---

## 📋 Manual Deployment Steps

### Step 1: Create Railway Project

1. Go to [Railway Dashboard](https://railway.app/dashboard)
2. Click **"New Project"**
3. Choose **"Deploy from GitHub repo"**
4. Select your `crios-dr-crowe-coder` repository
5. Railway will create an empty project

### Step 2: Add PostgreSQL Database

1. In your Railway project, click **"+ New"**
2. Select **"Database"** → **"Add PostgreSQL"**
3. Railway automatically provisions PostgreSQL 16
4. Note: `DATABASE_URL` is automatically set
5. Rename service to `Postgres` for clarity

**Configuration:**
- **Version**: PostgreSQL 16
- **Volume Size**: 10GB (adjust as needed)
- **Backups**: Enabled by default

### Step 3: Add Redis Cache

1. Click **"+ New"** again
2. Select **"Database"** → **"Add Redis"**
3. Railway automatically provisions Redis 7
4. Note: `REDIS_URL` is automatically set
5. Rename service to `Redis` for clarity

**Configuration:**
- **Version**: Redis 7
- **Volume Size**: 5GB
- **Persistence**: Enabled (AOF)

### Step 4: Deploy Backend Service

1. Click **"+ New"** → **"GitHub Repo"**
2. Select your repository
3. Click **"Add variables"** and configure:

```bash
# Reference PostgreSQL and Redis services
DATABASE_URL=${{Postgres.DATABASE_URL}}
REDIS_URL=${{Redis.REDIS_URL}}

# API Keys
ANTHROPIC_API_KEY=<your-anthropic-api-key>

# Application settings
ENV=production
PYTHONPATH=/app

# Security
ALLOWED_HOSTS=*
SECRET_KEY=<generate-with: openssl rand -base64 32>
```

4. Go to **"Settings"** → **"Build & Deploy"**:
   - **Root Directory**: Leave blank (use repository root)
   - **Dockerfile Path**: `platform/backend/Dockerfile`
   - **Build Command**: (leave default)

5. Go to **"Settings"** → **"Networking"**:
   - Enable **"Public Networking"**
   - Copy the public URL (e.g., `backend-production-xxxx.up.railway.app`)

6. Set **Health Check** (Settings → Health Check):
   - **Path**: `/health`
   - **Timeout**: 30 seconds

7. Click **"Deploy"** and wait for build to complete

### Step 5: Deploy Frontend Service

1. Click **"+ New"** → **"GitHub Repo"**
2. Select your repository again
3. Click **"Add variables"** and configure:

```bash
# Reference Backend service
NEXT_PUBLIC_API_URL=https://${{Backend.RAILWAY_PUBLIC_DOMAIN}}
NEXT_PUBLIC_WS_URL=wss://${{Backend.RAILWAY_PUBLIC_DOMAIN}}

# Next.js settings
NODE_ENV=production
```

4. Go to **"Settings"** → **"Build & Deploy"**:
   - **Root Directory**: `platform/frontend`
   - **Dockerfile Path**: `Dockerfile`

5. Go to **"Settings"** → **"Networking"**:
   - Enable **"Public Networking"**
   - Optionally add custom domain

6. Click **"Deploy"**

### Step 6: Deploy Celery Worker Service

1. Click **"+ New"** → **"GitHub Repo"**
2. Select your repository again
3. Click **"Add variables"**:

```bash
# Reference services
DATABASE_URL=${{Postgres.DATABASE_URL}}
REDIS_URL=${{Redis.REDIS_URL}}

# API Keys
ANTHROPIC_API_KEY=<your-anthropic-api-key>

# Celery settings
CELERY_BROKER_URL=${{Redis.REDIS_URL}}
CELERY_RESULT_BACKEND=${{Redis.REDIS_URL}}
C_FORCE_ROOT=true
PYTHONPATH=/app
```

4. Go to **"Settings"** → **"Build & Deploy"**:
   - **Root Directory**: Leave blank
   - **Dockerfile Path**: `platform/backend/Dockerfile`
   - **Custom Start Command**: `celery -A app.celery worker --loglevel=info --concurrency=2`

5. **Disable Public Networking** (this is a background worker)

6. Click **"Deploy"**

---

## 🔧 Post-Deployment Configuration

### Update Backend CORS Settings

After frontend deployment, update the Backend service environment variables:

```bash
CORS_ORIGINS=https://${{Frontend.RAILWAY_PUBLIC_DOMAIN}}
```

Redeploy the Backend service for changes to take effect.

### Verify Deployment

Test each service:

```bash
# Backend health check
curl https://your-backend-url.railway.app/health

# Frontend health check
curl https://your-frontend-url.railway.app/

# API Documentation
open https://your-backend-url.railway.app/docs
```

### Configure Custom Domains (Optional)

1. Go to Frontend service → **"Settings"** → **"Domains"**
2. Click **"Add Domain"**
3. Enter your custom domain (e.g., `crios.yourdomain.com`)
4. Add CNAME record in your DNS provider:
   - **Name**: `crios`
   - **Value**: `your-frontend-url.railway.app`

5. Repeat for Backend if needed (e.g., `api.yourdomain.com`)

---

## 📊 Service Architecture on Railway

```
┌─────────────────────────────────────────────────────────┐
│                    Railway Project                       │
│                                                          │
│  ┌──────────────┐         ┌──────────────┐             │
│  │  PostgreSQL  │◄────────┤   Backend    │             │
│  │  (Database)  │         │   (FastAPI)  │◄────┐       │
│  └──────────────┘         └──────┬───────┘     │       │
│                                   │             │       │
│  ┌──────────────┐                 │             │       │
│  │    Redis     │◄────────────────┤             │       │
│  │   (Cache)    │                 │             │       │
│  └──────────────┘                 │             │       │
│         ▲                         │             │       │
│         │                         │             │       │
│         │           ┌─────────────▼─────┐       │       │
│         └───────────┤ Celery Worker     │       │       │
│                     │ (Background Jobs) │       │       │
│                     └───────────────────┘       │       │
│                                                 │       │
│  ┌────────────────────────────────────┐        │       │
│  │          Frontend                  │────────┘       │
│  │         (Next.js)                  │                │
│  │  https://crios.railway.app         │                │
│  └────────────────────────────────────┘                │
│                                                         │
└─────────────────────────────────────────────────────────┘

Public Access:
  ├─ Frontend: https://your-app.railway.app
  └─ Backend API: https://your-backend.railway.app

Internal Communication:
  └─ Services use Railway's private networking
```

---

## 💰 Cost Estimation

Railway pricing (as of 2025):

| Service | Plan | Monthly Cost |
|---------|------|--------------|
| PostgreSQL 16 | 10GB volume | $5-10 |
| Redis 7 | 5GB volume | $3-5 |
| Backend | 1 vCPU, 1GB RAM | $5-10 |
| Frontend | 1 vCPU, 1GB RAM | $5-10 |
| Celery Worker | 1 vCPU, 1GB RAM | $5-10 |
| **Total** | **Estimated** | **$23-45/month** |

**Notes:**
- Railway includes $5 free credit per month
- Pricing based on usage (CPU, memory, network)
- First 100GB egress free per month
- Sleeping not recommended for production apps

---

## 🔒 Security Best Practices

### Environment Variables

1. **Never commit secrets** to git
2. Use Railway's environment variable management
3. Rotate secrets regularly:
   ```bash
   openssl rand -base64 32  # Generate new secret
   ```

### Database Security

```bash
# PostgreSQL: Enable SSL (automatic in Railway)
# Redis: Set password (automatic in Railway)
# Restrict network access to private only for databases
```

### API Security

Update `platform/backend/main.py` CORS settings:

```python
from fastapi.middleware.cors import CORSMiddleware

app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "https://your-frontend.railway.app",
        "https://your-custom-domain.com"
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
```

---

## 📈 Monitoring & Logging

### Railway Observability

1. **Logs**: Available in each service's **"Logs"** tab
2. **Metrics**: CPU, Memory, Network in **"Metrics"** tab
3. **Deployments**: History in **"Deployments"** tab

### Custom Monitoring

Add health checks to your backend:

```python
# platform/backend/app.py
from fastapi import FastAPI
from datetime import datetime

@app.get("/health")
async def health_check():
    return {
        "status": "healthy",
        "timestamp": datetime.utcnow().isoformat(),
        "service": "crios-backend"
    }

@app.get("/metrics")
async def metrics():
    # Add custom metrics
    return {
        "database": "connected",
        "redis": "connected",
        "version": "1.0.0"
    }
```

### External Monitoring (Optional)

- **Uptime Monitoring**: UptimeRobot, Pingdom
- **Error Tracking**: Sentry
- **APM**: New Relic, Datadog

---

## 🔄 CI/CD with Railway

Railway automatically deploys on git push to main branch.

### Custom Build Configuration

Create `railway.toml` in repository root (already included):

```toml
[build]
builder = "DOCKERFILE"
dockerfilePath = "Dockerfile"

[deploy]
numReplicas = 1
restartPolicyType = "ON_FAILURE"
restartPolicyMaxRetries = 10
```

### Deployment Workflow

```bash
# 1. Commit changes
git add .
git commit -m "feat: add new discovery algorithm"

# 2. Push to trigger deployment
git push origin main

# 3. Railway automatically:
#    - Detects changes
#    - Builds Docker image
#    - Runs health checks
#    - Deploys with zero-downtime
```

### Branch Deployments

Deploy feature branches:

1. Go to Railway service **"Settings"** → **"Source"**
2. Change **"Production Branch"** or add **"PR Deploys"**
3. Each PR gets a unique URL for testing

---

## 🛠️ Troubleshooting

### Common Issues

#### 1. Build Fails with "Cannot find Dockerfile"

**Solution:** Check Dockerfile path in Railway settings
```bash
# Backend: platform/backend/Dockerfile
# Frontend: Dockerfile (relative to root directory set to platform/frontend)
```

#### 2. Backend Can't Connect to Database

**Solution:** Verify environment variable references
```bash
# Correct syntax:
DATABASE_URL=${{Postgres.DATABASE_URL}}

# Not:
DATABASE_URL=${{POSTGRES.DATABASE_URL}}  # Wrong case
```

#### 3. Frontend Can't Reach Backend API

**Solution:** Check CORS and API URL
```bash
# Frontend env:
NEXT_PUBLIC_API_URL=https://${{Backend.RAILWAY_PUBLIC_DOMAIN}}

# Backend env (add):
CORS_ORIGINS=https://${{Frontend.RAILWAY_PUBLIC_DOMAIN}}
```

#### 4. Celery Worker Not Processing Tasks

**Solution:** Verify Redis connection and broker URL
```bash
# Check logs for connection errors
# Ensure CELERY_BROKER_URL points to Redis
CELERY_BROKER_URL=${{Redis.REDIS_URL}}
```

#### 5. High Memory Usage / OOM Kills

**Solution:** Optimize resources or upgrade
```bash
# Railway Settings → Resources
# Increase Memory limit (costs more)
# Or optimize Python/Node.js memory usage
```

### Debug Commands

```bash
# View service logs
railway logs --service backend

# SSH into container (Railway doesn't support, use logs instead)

# Check environment variables
railway variables --service backend

# Restart service
railway restart --service backend
```

---

## 📚 Additional Resources

- **Railway Documentation**: https://docs.railway.app
- **Railway Templates**: https://railway.app/templates
- **Railway Discord**: https://discord.gg/railway
- **CriOS Documentation**: See [README.md](README.md)
- **API Documentation**: Available at `https://your-backend.railway.app/docs`

---

## 🎉 Success Checklist

After deployment, verify:

- [ ] PostgreSQL database is running
- [ ] Redis cache is running
- [ ] Backend health check returns 200: `/health`
- [ ] Frontend loads successfully
- [ ] Celery worker is processing tasks
- [ ] API documentation accessible: `/docs`
- [ ] Molecular validation endpoint works: `POST /validate`
- [ ] Frontend can communicate with backend
- [ ] Environment variables properly set
- [ ] Custom domains configured (if applicable)
- [ ] SSL certificates active (automatic)
- [ ] Monitoring and logging working

---

## 🚦 Going to Production

### Before Production Launch

1. **Security Audit**
   - [ ] Review all environment variables
   - [ ] Enable rate limiting on API
   - [ ] Set up proper CORS origins
   - [ ] Use strong database passwords

2. **Performance Testing**
   - [ ] Load test API endpoints
   - [ ] Test Celery worker capacity
   - [ ] Verify database query performance
   - [ ] Check frontend load times

3. **Backup Strategy**
   - [ ] Configure PostgreSQL backups (Railway provides automatic backups)
   - [ ] Test restore procedures
   - [ ] Document recovery process

4. **Monitoring Setup**
   - [ ] Set up uptime monitoring
   - [ ] Configure error alerting
   - [ ] Set up log aggregation

### Production Environment Variables

See [`.env.railway.example`](.env.railway.example) for complete template.

---

## 📞 Support

- **Railway Support**: support@railway.app
- **CriOS Issues**: [GitHub Issues](https://github.com/michaelcrowe11/crios/issues)
- **Email**: michael@crowelogic.com

---

**CriOS Discovery Engine** • Deployed with Railway • Science before status. Discovery before profit.

© 2025 Crowe BioSystems
